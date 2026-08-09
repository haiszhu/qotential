/*
 * Worker-side progress decoding, formatting, and wall-clock throttling.
 *
 * The real Fortran core emits compact numeric events (stage code plus four i64
 * payloads and one f64) through the EM_JS bridge.  This module turns those into
 * make-style log lines and applies a ~250 ms wall-clock throttle to the
 * deterministic 2-percent milestones the Fortran side already produces.  It
 * computes no solver quantity; every displayed number originates in the solver.
 */
import type { SolverWasmModule } from './wasm-loader';

export type RawProgressCallback = NonNullable<SolverWasmModule['onSolverProgress']>;

export interface FormattedProgress {
  line: string;
  /** Stable per-phase key; begin and progress codes of one phase share it. */
  stageKey: string;
  /** Stage boundaries, 100-percent, scatter, render, and result bypass throttle. */
  alwaysForward: boolean;
}

const THROTTLE_MS = 250;

const integerFormat = new Intl.NumberFormat('en-US', { maximumFractionDigits: 0 });

/** Convert a solver count to a number only within the JS safe-integer range. */
export function safeCount(value: bigint, label: string): number {
  if (value < 0n || value > BigInt(Number.MAX_SAFE_INTEGER)) {
    throw new RangeError(`${label} count ${value} is outside the safe-integer range`);
  }
  return Number(value);
}

function fmt(value: bigint, label: string): string {
  return integerFormat.format(safeCount(value, label));
}

/** Percentage from exact bigint arithmetic, floored, matching the Fortran side. */
function percent(current: bigint, total: bigint): number {
  return Number((100n * current) / total);
}

function pad3(pct: number): string {
  return String(pct).padStart(3, ' ');
}

/** Scientific notation with a signed, zero-padded two-digit exponent. */
function formatGrfError(value: number): string {
  const text = value.toExponential(6);
  const [mantissa, rawExp] = text.split('e');
  const sign = rawExp.startsWith('-') ? '-' : '+';
  const digits = rawExp.replace(/^[+-]/, '').padStart(2, '0');
  return `${mantissa}e${sign}${digits}`;
}

/** A counted progress/100-percent line: `[label]  NN%  noun c / t`. */
function countedLine(
  label: string,
  noun: string,
  current: bigint,
  total: bigint,
): string {
  const pct = pad3(percent(current, total));
  return `[${label}] ${pct}%  ${noun} ${fmt(current, label)} / ${fmt(total, label)}`;
}

function progressStage(
  label: string,
  noun: string,
  current: bigint,
  total: bigint,
): FormattedProgress {
  return {
    line: countedLine(label, noun, current, total),
    stageKey: label,
    // The single 100-percent event must always reach the log.
    alwaysForward: current >= total,
  };
}

/**
 * Decode one raw progress event into a display line, or `undefined` for an
 * unknown stage code (which must be ignored, never fail the solve).  A count
 * outside the safe-integer range throws and is caught by the forwarder.
 */
export function formatProgressEvent(
  stage: number,
  current: bigint,
  total: bigint,
  aux0: bigint,
  _aux1: bigint,
  value: number,
): FormattedProgress | undefined {
  switch (stage) {
    case 1: {
      const isW7x = aux0 === 1n;
      const line = isW7x
        ? '[geometry] building W7-X adaptive charts; refined meshes may be quiet here for tens of seconds'
        : '[geometry] building the built-in surface mesh';
      return { line, stageKey: 'geometry', alwaysForward: true };
    }
    case 2:
      return {
        line: `[geometry] ${fmt(current, 'panels')} panels, ${fmt(total, 'source nodes')} source nodes`,
        stageKey: 'geometry',
        alwaysForward: true,
      };
    case 3:
      if (aux0 === 0n) {
        return {
          line: '[direct] evaluating naive Laplace GRF',
          stageKey: 'direct',
          alwaysForward: true,
        };
      }
      if (aux0 === 1n) {
        return {
          line: '[fmm] evaluating Laplace GRF with FMM3D',
          stageKey: 'fmm',
          alwaysForward: true,
        };
      }
      return undefined;
    case 4:
      if (aux0 === 0n) return progressStage('direct', 'block', current, total);
      if (aux0 === 1n) {
        if (!Number.isFinite(value) || value < 0) return undefined;
        return {
          line: `[fmm] completed Laplace GRF in ${value.toFixed(3)} s`,
          stageKey: 'fmm',
          alwaysForward: true,
        };
      }
      return undefined;
    case 5:
      return {
        line: '[close/count] finding close targets',
        stageKey: 'close/count',
        alwaysForward: true,
      };
    case 6:
      return progressStage('close/count', 'patch', current, total);
    case 7:
      return {
        line: '[close/rrq] applying close-panel correction',
        stageKey: 'close/rrq',
        alwaysForward: true,
      };
    case 8:
      return progressStage('close/rrq', 'patch', current, total);
    case 9:
      return {
        line: '[close/scatter] accumulating corrected values',
        stageKey: 'close/scatter',
        alwaysForward: true,
      };
    case 10:
      return {
        line: '[render] building visualization lattice',
        stageKey: 'render',
        alwaysForward: true,
      };
    case 11:
      return {
        line: `[result] GRF max relative error = ${formatGrfError(value)}`,
        stageKey: 'result',
        alwaysForward: true,
      };
    default:
      return undefined;
  }
}

export interface ProgressForwarderOptions {
  /** Receives each forwarded log line, in arrival order. */
  sink: (line: string) => void;
  /** Injectable monotonic clock in milliseconds; defaults to performance.now. */
  now?: () => number;
}

/**
 * Build the raw callback installed on `Module.onSolverProgress`.  It decodes,
 * applies the per-stage 250 ms throttle to ordinary milestones, and forwards
 * surviving lines to `sink`.  It never throws: a decode or sink failure is
 * swallowed so progress reporting cannot perturb the solve.
 */
export function createProgressForwarder(
  options: ProgressForwarderOptions,
): RawProgressCallback {
  const now = options.now ?? (() => performance.now());
  const lastForwarded = new Map<string, number>();

  return (stage, current, total, aux0, aux1, value) => {
    try {
      const event = formatProgressEvent(stage, current, total, aux0, aux1, value);
      if (!event) return;

      const at = now();
      if (!event.alwaysForward) {
        const previous = lastForwarded.get(event.stageKey);
        if (previous !== undefined && at - previous < THROTTLE_MS) return;
      }
      lastForwarded.set(event.stageKey, at);
      options.sink(event.line);
    } catch {
      // Observational only: never let progress handling abort the solve.
    }
  };
}

/**
 * Install `callback` on the module for exactly the duration of `run`, then
 * remove it — after both success and failure.  The module stays reusable for a
 * later sequential solve.
 */
export function withSolverProgress<T>(
  wasm: SolverWasmModule,
  callback: RawProgressCallback,
  run: () => T,
): T {
  wasm.onSolverProgress = callback;
  try {
    return run();
  } finally {
    delete wasm.onSolverProgress;
  }
}

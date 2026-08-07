export const LOG_ERROR_DOMAIN = [-16, 0] as const;

const STOPS = [
  [0.035, 0.055, 0.18],
  [0.02, 0.42, 0.67],
  [0.10, 0.78, 0.69],
  [0.95, 0.79, 0.24],
  [0.78, 0.12, 0.14],
] as const;

export function colorizeLogError(values: Float32Array): Float32Array {
  const colors = new Float32Array(4 * values.length);
  const [lo, hi] = LOG_ERROR_DOMAIN;
  for (let i = 0; i < values.length; ++i) {
    if (!Number.isFinite(values[i])) throw new Error(`log error value ${i} is not finite`);
    const unit = Math.max(0, Math.min(1, (values[i] - lo) / (hi - lo)));
    const scaled = unit * (STOPS.length - 1);
    const left = Math.min(STOPS.length - 2, Math.floor(scaled));
    const mix = scaled - left;
    for (let channel = 0; channel < 3; ++channel) {
      colors[4 * i + channel] = STOPS[left][channel] * (1 - mix) +
        STOPS[left + 1][channel] * mix;
    }
    colors[4 * i + 3] = 1;
  }
  return colors;
}

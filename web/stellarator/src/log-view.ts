/*
 * Pure model for the append-only run log.
 *
 * Only string state and the scroll-follow decision live here so they can be
 * unit-tested without a DOM.  The browser layer (main.ts) owns the actual log
 * element and applies these decisions to it.
 */

export interface LogBuffer {
  clear(): void;
  append(line: string): void;
  text(): string;
  lines(): readonly string[];
}

export function createLogBuffer(): LogBuffer {
  let lines: string[] = [];
  return {
    clear(): void {
      lines = [];
    },
    append(line: string): void {
      lines.push(line);
    },
    text(): string {
      return lines.join('\n');
    },
    lines(): readonly string[] {
      return lines;
    },
  };
}

/**
 * Whether the log should auto-scroll to the newest line.  True when the
 * viewport is already within `tolerance` pixels of the bottom (or the content
 * is short enough to be fully visible); false once the user has scrolled up.
 */
export function shouldFollow(
  scrollTop: number,
  clientHeight: number,
  scrollHeight: number,
  tolerance = 8,
): boolean {
  return scrollHeight - (scrollTop + clientHeight) <= tolerance;
}

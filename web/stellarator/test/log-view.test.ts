import { describe, expect, it } from 'vitest';
import { createLogBuffer, shouldFollow } from '../src/log-view';

describe('createLogBuffer', () => {
  it('preserves arrival order of appended messages', () => {
    const buffer = createLogBuffer();
    buffer.append('[geometry] building');
    buffer.append('[direct] 100%  block 1,805 / 1,805');
    buffer.append('[result] completed in 1.0 s');
    expect(buffer.lines()).toEqual([
      '[geometry] building',
      '[direct] 100%  block 1,805 / 1,805',
      '[result] completed in 1.0 s',
    ]);
    expect(buffer.text()).toBe(
      '[geometry] building\n[direct] 100%  block 1,805 / 1,805\n[result] completed in 1.0 s',
    );
  });

  it('clears lines from the previous run', () => {
    const buffer = createLogBuffer();
    buffer.append('[cancelled] solve terminated by user');
    buffer.clear();
    buffer.append('[geometry] building the built-in surface mesh');
    expect(buffer.lines()).toEqual(['[geometry] building the built-in surface mesh']);
  });

  it('keeps completion, error, and cancellation lines until the next clear', () => {
    const buffer = createLogBuffer();
    buffer.append('[error] status=215, solver error=215');
    expect(buffer.text()).toBe('[error] status=215, solver error=215');
    // Still present without a clear.
    buffer.append('[cancelled] solve terminated by user');
    expect(buffer.lines()).toHaveLength(2);
  });
});

describe('shouldFollow', () => {
  it('follows when the viewport is within tolerance of the bottom', () => {
    // scrollHeight 1000, viewport 200 tall, scrolled to 792 -> 8px from bottom.
    expect(shouldFollow(792, 200, 1000)).toBe(true);
  });

  it('does not follow when the user has scrolled upward', () => {
    expect(shouldFollow(400, 200, 1000)).toBe(false);
  });

  it('follows an empty or short log fully within the viewport', () => {
    expect(shouldFollow(0, 200, 0)).toBe(true);
    expect(shouldFollow(0, 200, 150)).toBe(true);
  });

  it('respects a custom tolerance', () => {
    expect(shouldFollow(700, 200, 1000, 8)).toBe(false);
    expect(shouldFollow(700, 200, 1000, 100)).toBe(true);
  });
});

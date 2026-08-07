import { colorizeLogError } from './colormap';

export interface PreparedGpuMesh {
  positions: Float32Array;
  colors: Float32Array;
  triangles: Uint32Array;
}

export function prepareGpuMesh(
  positions: Float32Array,
  logError: Float32Array,
  triangles: Uint32Array,
): PreparedGpuMesh {
  if (positions.length !== 3 * logError.length) throw new Error('mesh shape mismatch');
  const normalized = positions.slice();
  const lo = [Infinity, Infinity, Infinity];
  const hi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < logError.length; ++i) {
    for (let axis = 0; axis < 3; ++axis) {
      const value = positions[3 * i + axis];
      lo[axis] = Math.min(lo[axis], value);
      hi[axis] = Math.max(hi[axis], value);
    }
  }
  const center = lo.map((value, axis) => 0.5 * (value + hi[axis]));
  let radius = 0;
  for (let i = 0; i < logError.length; ++i) {
    for (let axis = 0; axis < 3; ++axis) {
      radius = Math.max(radius, Math.abs(positions[3 * i + axis] - center[axis]));
    }
  }
  if (!(radius > 0)) throw new Error('mesh has zero spatial extent');
  for (let i = 0; i < normalized.length; ++i) {
    normalized[i] = (normalized[i] - center[i % 3]) / radius;
  }
  return {
    positions: normalized,
    colors: colorizeLogError(logError),
    triangles: triangles.slice(),
  };
}

const SHADER = `
struct Uniforms { mvp: mat4x4<f32> };
@group(0) @binding(0) var<uniform> uniforms: Uniforms;

struct VertexOut {
  @builtin(position) position: vec4<f32>,
  @location(0) color: vec4<f32>,
};

@vertex fn vertex_main(
  @location(0) position: vec3<f32>,
  @location(1) color: vec4<f32>,
) -> VertexOut {
  var out: VertexOut;
  out.position = uniforms.mvp * vec4<f32>(position, 1.0);
  out.color = color;
  return out;
}

@fragment fn fragment_main(in: VertexOut) -> @location(0) vec4<f32> {
  return in.color;
}
`;

function multiply(a: Float32Array, b: Float32Array): Float32Array {
  const out = new Float32Array(16);
  for (let col = 0; col < 4; ++col) {
    for (let row = 0; row < 4; ++row) {
      for (let k = 0; k < 4; ++k) out[col * 4 + row] += a[k * 4 + row] * b[col * 4 + k];
    }
  }
  return out;
}

function cameraMatrix(yaw: number, pitch: number, zoom: number, aspect: number): Float32Array {
  const f = 1 / Math.tan(0.5 * Math.PI / 4);
  const projection = new Float32Array([
    f / aspect, 0, 0, 0,
    0, f, 0, 0,
    0, 0, -1.002, -1,
    0, 0, -0.2002, 0,
  ]);
  const cy = Math.cos(yaw), sy = Math.sin(yaw);
  const cp = Math.cos(pitch), sp = Math.sin(pitch);
  const rotationY = new Float32Array([
    cy, 0, -sy, 0, 0, 1, 0, 0, sy, 0, cy, 0, 0, 0, 0, 1,
  ]);
  const rotationX = new Float32Array([
    1, 0, 0, 0, 0, cp, sp, 0, 0, -sp, cp, 0, 0, 0, 0, 1,
  ]);
  const view = new Float32Array([
    1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, -zoom, 1,
  ]);
  return multiply(projection, multiply(view, multiply(rotationX, rotationY)));
}

export class StellaratorRenderer {
  private readonly context: GPUCanvasContext;
  private readonly format: GPUTextureFormat;
  private readonly pipeline: GPURenderPipeline;
  private readonly uniformBuffer: GPUBuffer;
  private readonly bindGroup: GPUBindGroup;
  private positionBuffer?: GPUBuffer;
  private colorBuffer?: GPUBuffer;
  private indexBuffer?: GPUBuffer;
  private depthTexture?: GPUTexture;
  private indexCount = 0;
  private yaw = -0.6;
  private pitch = 0.45;
  private zoom = 3.3;

  private constructor(
    private readonly canvas: HTMLCanvasElement,
    private readonly device: GPUDevice,
  ) {
    const context = canvas.getContext('webgpu');
    if (!context) throw new Error('WebGPU canvas context is unavailable');
    this.context = context;
    this.format = navigator.gpu.getPreferredCanvasFormat();
    context.configure({ device, format: this.format, alphaMode: 'premultiplied' });
    const shader = device.createShaderModule({ code: SHADER });
    this.pipeline = device.createRenderPipeline({
      layout: 'auto',
      vertex: {
        module: shader,
        entryPoint: 'vertex_main',
        buffers: [
          { arrayStride: 12, attributes: [{ shaderLocation: 0, offset: 0, format: 'float32x3' }] },
          { arrayStride: 16, attributes: [{ shaderLocation: 1, offset: 0, format: 'float32x4' }] },
        ],
      },
      fragment: { module: shader, entryPoint: 'fragment_main', targets: [{ format: this.format }] },
      primitive: { topology: 'triangle-list', cullMode: 'none' },
      depthStencil: { format: 'depth24plus', depthWriteEnabled: true, depthCompare: 'less' },
    });
    this.uniformBuffer = device.createBuffer({
      size: 64,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
    });
    this.bindGroup = device.createBindGroup({
      layout: this.pipeline.getBindGroupLayout(0),
      entries: [{ binding: 0, resource: { buffer: this.uniformBuffer } }],
    });
    this.installControls();
    new ResizeObserver(() => this.draw()).observe(canvas);
  }

  static async create(canvas: HTMLCanvasElement): Promise<StellaratorRenderer> {
    if (!navigator.gpu) throw new Error('WebGPU is not available in this browser');
    const adapter = await navigator.gpu.requestAdapter();
    if (!adapter) throw new Error('No WebGPU adapter is available');
    return new StellaratorRenderer(canvas, await adapter.requestDevice());
  }

  setMesh(mesh: PreparedGpuMesh): void {
    this.positionBuffer?.destroy();
    this.colorBuffer?.destroy();
    this.indexBuffer?.destroy();
    this.positionBuffer = this.buffer(mesh.positions, GPUBufferUsage.VERTEX);
    this.colorBuffer = this.buffer(mesh.colors, GPUBufferUsage.VERTEX);
    this.indexBuffer = this.buffer(mesh.triangles, GPUBufferUsage.INDEX);
    this.indexCount = mesh.triangles.length;
    this.draw();
  }

  private buffer(data: Float32Array | Uint32Array, usage: GPUBufferUsageFlags): GPUBuffer {
    const buffer = this.device.createBuffer({
      size: Math.max(4, (data.byteLength + 3) & ~3),
      usage: usage | GPUBufferUsage.COPY_DST,
    });
    this.device.queue.writeBuffer(
      buffer,
      0,
      data.buffer as ArrayBuffer,
      data.byteOffset,
      data.byteLength,
    );
    return buffer;
  }

  private installControls(): void {
    let dragging = false, x = 0, y = 0;
    this.canvas.addEventListener('pointerdown', (event) => {
      dragging = true; x = event.clientX; y = event.clientY;
      this.canvas.setPointerCapture(event.pointerId);
    });
    this.canvas.addEventListener('pointermove', (event) => {
      if (!dragging) return;
      this.yaw += (event.clientX - x) * 0.008;
      this.pitch = Math.max(-1.35, Math.min(1.35, this.pitch + (event.clientY - y) * 0.008));
      x = event.clientX; y = event.clientY; this.draw();
    });
    this.canvas.addEventListener('pointerup', () => { dragging = false; });
    this.canvas.addEventListener('wheel', (event) => {
      event.preventDefault();
      this.zoom = Math.max(2.1, Math.min(7, this.zoom * Math.exp(event.deltaY * 0.001)));
      this.draw();
    }, { passive: false });
  }

  draw(): void {
    if (!this.indexBuffer || !this.positionBuffer || !this.colorBuffer || !this.indexCount) return;
    const ratio = Math.min(2, window.devicePixelRatio || 1);
    const width = Math.max(1, Math.floor(this.canvas.clientWidth * ratio));
    const height = Math.max(1, Math.floor(this.canvas.clientHeight * ratio));
    if (!this.depthTexture || this.canvas.width !== width || this.canvas.height !== height) {
      this.canvas.width = width; this.canvas.height = height;
      this.depthTexture?.destroy();
      this.depthTexture = this.device.createTexture({
        size: [width, height], format: 'depth24plus', usage: GPUTextureUsage.RENDER_ATTACHMENT,
      });
    }
    const camera = cameraMatrix(this.yaw, this.pitch, this.zoom, width / height);
    this.device.queue.writeBuffer(
      this.uniformBuffer,
      0,
      camera.buffer as ArrayBuffer,
      camera.byteOffset,
      camera.byteLength,
    );
    const encoder = this.device.createCommandEncoder();
    const pass = encoder.beginRenderPass({
      colorAttachments: [{
        view: this.context.getCurrentTexture().createView(),
        clearValue: { r: 0.012, g: 0.018, b: 0.035, a: 1 },
        loadOp: 'clear', storeOp: 'store',
      }],
      depthStencilAttachment: {
        view: this.depthTexture!.createView(),
        depthClearValue: 1, depthLoadOp: 'clear', depthStoreOp: 'store',
      },
    });
    pass.setPipeline(this.pipeline);
    pass.setBindGroup(0, this.bindGroup);
    pass.setVertexBuffer(0, this.positionBuffer);
    pass.setVertexBuffer(1, this.colorBuffer);
    pass.setIndexBuffer(this.indexBuffer, 'uint32');
    pass.drawIndexed(this.indexCount);
    pass.end();
    this.device.queue.submit([encoder.finish()]);
  }
}

// Standalone CPU Wasm: explicit imports, native exports, and fresh heap views.
// Browser/GPU control stays in JS; this module has no emulated filesystem.
export async function createCpuWasm(source, {imports = {}, print = console.log, printErr = console.error} = {}) {
  const compiled = source instanceof WebAssembly.Module ? source :
    await WebAssembly.compile(source instanceof Uint8Array || source instanceof ArrayBuffer ? source :
      await (async () => {
        const response = await fetch(source);
        if (!response.ok) throw Error(`Could not load CPU Wasm: ${response.status}`);
        return response.arrayBuffer();
      })());
  let memory;
  const writeZero = pointer => new DataView(memory.buffer).setUint32(pointer >>> 0, 0, true);
  const messages = new Map(), decoders = new Map();
  const host = {
    env: {emscripten_notify_memory_growth() {}},
    wasi_snapshot_preview1: {
      fd_close: () => 8,
      clock_time_get(clock, precision, pointer) {
        const now = clock === 0 ? Date.now() : performance.now();
        new DataView(memory.buffer).setBigUint64(pointer >>> 0, BigInt(Math.floor(now * 1e6)), true);
        return 0;
      },
      fd_read(fd, vectors, count, written) { writeZero(written); return fd === 0 ? 0 : 8; },
      fd_write(fd, vectors, count, written) {
        if (fd !== 1 && fd !== 2) return 8;
        const heap = new Uint8Array(memory.buffer), words = new DataView(memory.buffer);
        if (!decoders.has(fd)) decoders.set(fd, new TextDecoder());
        let bytes = 0, line = messages.get(fd) || '';
        for (let i = 0; i < count; ++i) {
          const start = words.getUint32(vectors + 8*i, true), size = words.getUint32(vectors + 8*i + 4, true);
          line += decoders.get(fd).decode(heap.subarray(start, start + size), {stream:true}); bytes += size;
        }
        const lines = line.split('\n'); messages.set(fd, lines.pop());
        for (const value of lines) (fd === 1 ? print : printErr)(value);
        words.setUint32(written, bytes, true); return 0;
      },
      random_get(pointer, count) {
        const data = new Uint8Array(memory.buffer, pointer >>> 0, count);
        for (let i = 0; i < count; i += 65536) crypto.getRandomValues(data.subarray(i, i + 65536));
        return 0;
      },
      environ_sizes_get(count, size) { writeZero(count); writeZero(size); return 0; },
      environ_get(environment) { writeZero(environment); return 0; }
    }
  };
  for (const [name, values] of Object.entries(imports)) host[name] = {...host[name], ...values};
  for (const entry of WebAssembly.Module.imports(compiled))
    if (host[entry.module]?.[entry.name] === undefined)
      throw Error(`Unsupported CPU Wasm import: ${entry.module}.${entry.name}`);
  const instance = await WebAssembly.instantiate(compiled, host), exports = instance.exports;
  memory = exports.memory;
  if (!(memory instanceof WebAssembly.Memory) || typeof exports._initialize !== 'function')
    throw Error('CPU Wasm must export memory and _initialize.');
  exports._initialize();
  return cpuMemoryAccess(exports,memory);
}

// Both the standalone and pthread loaders expose the same CPU ABI.
export function cpuMemoryAccess(exports, memory) {
  const decoder = new TextDecoder(), encoder = new TextEncoder();
  return {
    exports, memory,
    view(Type, pointer = 0, count) { return new Type(memory.buffer, pointer >>> 0, count); },
    string(pointer) {
      const bytes = new Uint8Array(memory.buffer), start = pointer >>> 0, end = bytes.indexOf(0, start);
      if (!start || end < 0) throw Error('Invalid CPU Wasm string.');
      // Some browsers still reject shared views in TextDecoder. Copy only the
      // string; numerical tables and uniforms keep their borrowed heap views.
      return decoder.decode(bytes.buffer instanceof ArrayBuffer ? bytes.subarray(start,end) : bytes.slice(start,end));
    },
    withString(value, callback) {
      const bytes = encoder.encode(value), pointer = exports.malloc(bytes.length + 1);
      if (!pointer) throw Error('CPU Wasm allocation failed.');
      try {
        const heap = new Uint8Array(memory.buffer);
        heap.set(bytes, pointer >>> 0); heap[(pointer >>> 0) + bytes.length] = 0;
        const result = callback(pointer);
        if (result && typeof result.then === 'function') throw Error('withString requires a synchronous callback.');
        return result;
      } finally { exports.free(pointer); }
    }
  };
}

// File operations execute inside native WasmFS; JS only copies owned bytes.
export function cpuFiles(cpu) {
  const api = cpu.exports;
  const fail = () => Error(cpu.string(api.equilibrium_file_error()));
  return {
    write(path, bytes) {
      if (!(bytes instanceof Uint8Array) || bytes.length > 0x7fffffff)
        throw Error('Invalid CPU file bytes.');
      const pointer = api.malloc(Math.max(1, bytes.length));
      if (!pointer) throw Error('CPU Wasm allocation failed.');
      try {
        cpu.view(Uint8Array, pointer, bytes.length).set(bytes);
        cpu.withString(path, name => {
          if (api.equilibrium_file_write(name, pointer, bytes.length)) throw fail();
        });
      } finally { api.free(pointer); }
    },
    read(path) {
      return cpu.withString(path, name => {
        const pointer = api.equilibrium_file_read(name);
        if (!pointer) throw fail();
        return cpu.view(Uint8Array, pointer, api.equilibrium_file_size()).slice();
      });
    },
    remove(path) { cpu.withString(path, name => api.equilibrium_file_remove(name)); }
  };
}

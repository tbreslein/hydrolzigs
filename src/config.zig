pub const Config = struct {
    mesh: MeshConfig,

    pub fn validate(comptime self: Config) void {
        self.mesh.validate();
    }
};

const MeshConfig = struct {
    n_comp: u32,
    n_gc: u32 = 2,
    xi_in: f64,
    xi_out: f64,

    fn validate(comptime self: MeshConfig) void {
        if (self.xi_in >= self.xi_out) {
            @compileError("MeshConfig.xi_in < MeshConfig.xi_out must hold!");
        }
    }
};

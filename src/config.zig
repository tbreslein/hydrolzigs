pub const Config = struct {
    mesh: MeshConfig,

    pub fn validate(comptime self: Config) void {
        self.mesh.validate();
    }
};

pub const MeshConfig = struct {
    pub const Type = enum {
        cartesian,
    };

    type: Type,
    n: u32,
    xi_in: f64,
    xi_out: f64,

    fn validate(comptime self: MeshConfig) void {
        if (!(self.n > 5)) {
            @compileError("MeshConfig.n_all > 5 must hold!");
        }
        if (!(self.xi_in < self.xi_out)) {
            @compileError("MeshConfig.xi_in < MeshConfig.xi_out must hold!");
        }
    }
};

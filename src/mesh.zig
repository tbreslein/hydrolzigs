const Config = @import("config.zig").Config;

pub fn Mesh(comptime N: u32) type {
    return struct {
        comptime n_all: u32 = N,
        comptime n_comp: u32 = N - 2*2,
        comptime n_gc: u32 = 2,
        xi_in: f64,
        xi_out: f64,
    };
}

pub fn init_mesh(comptime conf: Config) Mesh(conf.mesh.n_comp + 2*2) {
    return Mesh(conf.mesh.n_comp + 2*2){
        .xi_in = 
    };
}

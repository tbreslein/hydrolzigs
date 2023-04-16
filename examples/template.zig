const std = @import("std");
const hydrol = @import("hydrol");

const conf = hydrol.Config{
    .mesh = .{
        .type = .cartesian,
        .n = 10,
        .xi_in = 1.0,
        .xi_out = 2.0,
    },
    .physics = .{
        .type = .euler1d_isothermal,
    },
};

pub fn main() !void {
    const mesh = hydrol.Mesh(conf){};
    var u = hydrol.Physics(conf){};
    hydrol.run(conf, &u, mesh);
}

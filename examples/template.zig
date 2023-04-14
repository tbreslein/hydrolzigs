const std = @import("std");
const hydrol = @import("hydrol");

const conf = hydrol.config.Config{
    .mesh = .{
        .n_comp = 100,
        .xi_in = 1.0,
        .xi_out = 2.0,
    },
};

pub fn main() !void {
    const n_gc = hydrol.init(conf);
    std.debug.print("n_gc = {}\n", .{n_gc});
    std.debug.print("conf = {}\n", .{conf});
}

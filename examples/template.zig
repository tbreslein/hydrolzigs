const std = @import("std");
const hydrol = @import("hydrol");

const conf = hydrol.config.Config{
    .mesh = .{
        .type = .cartesian,
        .n = 10,
        .xi_in = 1.0,
        .xi_out = 2.0,
    },
};

pub fn main() !void {
    const mesh = hydrol.init(conf);
    std.debug.print("mesh = {}\n", .{mesh});
    std.debug.print("conf = {}\n", .{conf});
}

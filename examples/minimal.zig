const std = @import("std");
const hydrol = @import("hydrol");

pub fn main() !void {
    // Prints to stderr (it's a shortcut based on `std.io.getStdErr()`)
    std.debug.print("All your {s} are belong to us!\n", .{"codebase"});
    const x = hydrol.add(2, 3);
    std.debug.print("add(2,3) = {}!\n", .{x});
}

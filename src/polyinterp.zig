fn polyInterp(allocator: mem.Allocator, x: f64, xs: []f64, ys: []f64, n: usize) !f64 {
    const offset: f64 = blk: {
        var min: f64 = math.inf(f64);
        var tmp: f64 = undefined;
        var ind: usize = undefined;
        for (0..n) |i| {
            tmp = @abs(x - xs[i]);
            if (tmp == 0.0) return ys[i];
            if (tmp < min) {
                min = tmp;
                ind = i;
            }
        }
        break :blk ys[ind];
    };

    const table: []f64 = try allocator.alloc(f64, n);
    defer allocator.free(table);

    for (0..n) |i| table[i] = ys[i] - offset;

    for (1..n) |j| {
        for (0..n - j, 1..) |i, ip1| {
            table[i] += (table[i] - table[ip1]) * (x - xs[i]) / (xs[i] - xs[i + j]);
        }
    }

    return offset + table[0];
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

fn func(x: f64) f64 {
    return (x + 1) * (x - 1) * (x - 2);
}

test "polyInterp" {
    const page = testing.allocator;
    var xs: [4]f64 = .{ -1.5, 0.0, 1.5, 2.5 };
    var ys: [4]f64 = .{ -4.375, 2.0, -0.625, 2.625 };
    try testing.expectEqual(func(-2.5), try polyInterp(page, -2.5, &xs, &ys, 4));
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

const std = @import("std");
const mem = std.mem;
const math = std.math;
const testing = std.testing;

const Array = @import("./array.zig").Array;

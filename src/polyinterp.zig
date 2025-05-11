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

test "polyInterp" {
    const page = testing.allocator;
    var xs: [4]f64 = .{ -1.5, 0.5, 1.5, 2.5 };
    var ys: [4]f64 = undefined;
    for (xs, &ys) |x, *y| y.* = func(x);
    try testing.expectEqual(func(-2.5), try polyInterp(page, -2.5, &xs, &ys, 4));
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

const NewtonPolynomial = struct {
    coeff: []f64,
    table: []f64,

    const Self: type = @This();

    fn init(allocator: mem.Allocator, order: usize) !*Self {
        const self: *Self = try allocator.create(Self);
        errdefer allocator.destroy(self);

        self.coeff = try allocator.alloc(f64, order);
        errdefer allocator.free(self.coeff);

        self.table = try allocator.alloc(f64, order);

        return self;
    }

    fn deinit(self: *const Self, allocator: mem.Allocator) void {
        allocator.free(self.table);
        allocator.free(self.coeff);
        allocator.destroy(self);
    }

    fn makepoly(self: *const Self, xarr: []f64, yarr: []f64) void {
        debug.assert(self.coeff.len == xarr.len);
        debug.assert(self.coeff.len == yarr.len);

        self.table[0] = yarr[0];
        self.coeff[0] = self.table[0];

        for (1..self.table.len) |k| {
            self.table[k] = yarr[k];

            var i: usize = k - 1;
            while (true) : (i -= 1) {
                self.table[i] = (self.table[i] - self.table[i + 1]) / (xarr[i] - xarr[k]);
                if (i == 0) break;
            }

            self.coeff[k] = self.table[0];
        }
    }

    fn evalpoly(self: *const Self, xarr: []f64, newx: f64) f64 {
        debug.assert(self.coeff.len == xarr.len);

        var n: usize = xarr.len - 1;
        var p: f64 = self.coeff[n];

        n -= 1;
        while (true) : (n -= 1) {
            p = p * (newx - xarr[n]) + self.coeff[n];
            if (n == 0) break;
        }

        return p;
    }
};

test "NewtonPolynomial" {
    const page = testing.allocator;

    var xs: [4]f64 = .{ -1.5, 0.5, 1.5, 2.5 };
    var ys: [4]f64 = undefined;
    for (xs, &ys) |x, *y| y.* = func(x);

    const poly: *NewtonPolynomial = try NewtonPolynomial.init(page, 4);
    defer poly.deinit(page);

    poly.makepoly(&xs, &ys);
    try testing.expectEqual(func(-2.5), poly.evalpoly(&xs, -2.5));
}

// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

fn func(x: f64) f64 {
    return (x + 1) * (x - 1) * (x - 2);
}

const std = @import("std");
const mem = std.mem;
const math = std.math;
const debug = std.debug;
const testing = std.testing;

const Array = @import("./array.zig").Array;

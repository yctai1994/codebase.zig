inline fn inv(x: f64) f64 {
    return 1.0 / x;
}

inline fn sqr(x: f64) f64 {
    return x * x;
}

inline fn pow(base: f64, expo: f64) f64 {
    return @exp(expo * @log(base));
}

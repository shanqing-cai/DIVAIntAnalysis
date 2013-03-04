function gen_rosenberg_source
F0 = 160;
N1 = 0.5;
N2 = 0.6;
fs = 16e3;
amp = 0.5;

NCycles= 1000;

wavfn = 'rosen_source.wav';

x = rosenberg(N1, N2, F0, fs);
w = repmat(x, 1, NCycles);
w = demean(amp * w);

wavwrite(w, fs, wavfn);

return
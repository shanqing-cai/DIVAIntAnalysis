function gen_pulse_source
F0 = 160;
fs = 16e3;
amp = 0.5;
dur = 3;

w = zeros(dur * fs, 1);
nspace = fs / F0;

w(1 : nspace : end) = amp;


wavfn = 'pulse_source.16k.wav';
wavwrite(w, fs, wavfn);

return
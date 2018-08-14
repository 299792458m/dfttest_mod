# dfttest
this is fork of dfttest 1.9.4. and some speed tune for my enviroment is added.
mainly tuned at dither function(it costs wastefully a lot).
And not tested a lot.

In my case(i7-4790K 4core8threads-CPU), the cost of dfttest.dll itself(not includeed fftw)
reduced 30% at dither=1(default) and dither>1 may also get good performance.
But dither=0 is little affected by this fork(because other tune is not critical).

I wish this fork become your help.

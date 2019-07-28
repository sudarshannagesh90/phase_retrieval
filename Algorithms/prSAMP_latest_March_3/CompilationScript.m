%mex -g -largeArrayDims -ldl CFLAGS="\$CFLAGS -std=c99" -output bin/swamp src/swamp.c src/channels/gaussian.c src/common/sort.c src/solvers/amp_alt.c src/solvers/gamp.c src/channels/nophase.c src/priors/binary.c src/solvers/amp.c src/channels/pm1.c src/priors/gb.c src/solvers/amp_dense.c 
% eval('mex -g -largeArrayDims -ldl CFLAGS="\$CFLAGS -std=c99" -output myout ./cswgamp.c ./solvers/cgamp.c ./gsl/* ./channels/* ./common/*')
% eval('mex -g -largeArrayDims -ldl CFLAGS="\$CFLAGS -std=c99" -output bin/myout ./channels/cpr.c ./common/sort.c ./gsl/bessel_I0.c ./gsl/bessel_I1.c ./gsl/cheb_eval.c ./priors/cgb.c ./solvers/cgamp.c ./cswgamp.c -lgsl -lgslcblas')
% eval('mex -g -largeArrayDims -ldl CFLAGS="\$CFLAGS -std=c99" -output bin/myprSwAMP ./channels/cpr.c ./common/sort.c ./gsl/bessel_I0.c ./gsl/bessel_I1.c ./priors/cgb.c ./solvers/cgamp.c ./cswgamp.c -lgsl -lgslcblas')

%The following works on Linux
%eval('mex -g -largeArrayDims -ldl CFLAGS="\$CFLAGS -std=c99" -output bin/myprSwAMP ./channels/cpr.c ./common/sort.c ./gsl/bessel_I0.c ./gsl/bessel_I1.c ./priors/cgb.c ./solvers/cgamp.c ./cswgamp.c')
%eval('mex -g -largeArrayDims -ldl CFLAGS="\$CFLAGS -std=c99" -output bin/PartialprSwAMP ./channels/cpr.c ./common/sort.c ./gsl/bessel_I0.c ./gsl/bessel_I1.c ./priors/cgb.c ./priors/svt.c ./solvers/partial_cgamp.c ./partial_cswgamp.c')

%eval('mex -g -largeArrayDims -ldl CFLAGS="\$CFLAGS -std=c99" -output bin/ArbMatprSwAMP ./channels/cpr.c ./common/sort.c ./gsl/bessel_I0.c ./gsl/bessel_I1.c ./priors/cgb.c ./priors/cgb_nzmean.c ./priors/cgb_nzmean_i.c ./solvers/arbmat_cgamp.c ./arbmat_cswgamp.c')

eval('mex -g -largeArrayDims -ldl CFLAGS="\$CFLAGS -std=c99" -output bin/prSAMP_newchannel ./channels/cpr.c ./common/sort.c ./gsl/bessel_I0.c ./gsl/bessel_I1.c ./priors/cgb.c ./priors/cgb_nzmean.c ./priors/cgb_nzmean_i.c ./solvers/arbmat_cgamp_newchannel.c ./arbmat_cswgamp_newchannel.c')
% warning('I never compiled cheb_eval.  Maybe it is not necessary');

% %The following should work, but doesn't, on Windows
% eval('mex CFLAGS="\$CFLAGS -std=c99" -output bin/myprSwAMP ./channels/cpr.c ./common/sort.c ./gsl/bessel_I0.c ./gsl/bessel_I1.c ./priors/cgb.c ./solvers/cgamp.c ./cswgamp.c')
% warning('I never compiled cheb_eval.  Maybe it is not necessary');


clear mex;

%%%--WMB FUNCTIONS--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mex -I../../cpp/include mexWmb.cpp -DIMPSAMP -output ../wmbImpSample
%mex -I../../cpp/include mexFactorIdx.cpp -DVALUE -output ../value      
%mex -I../../cpp/include mexFactorIdx.cpp -DSUBV2IND -output ../subv2ind


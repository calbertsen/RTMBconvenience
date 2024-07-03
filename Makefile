R?=R-4-2-2

roxygen:
	echo 'roxygen2::roxygenize("RTMBconvenience")' | $(R) --slave
	sed -i '/RcppExport void R_init/ s/^/void rtmb_set_shared_pointers();\n/' RTMBconvenience/src/RcppExports.cpp
	sed -i '/R_useDynamicSymbols/ s/$$/\n    rtmb_set_shared_pointers();/' RTMBconvenience/src/RcppExports.cpp

rcpp:
	echo 'Rcpp::compileAttributes("RTMBconvenience", verbose=TRUE)' | $(R) --slave
	sed -i '/RcppExport void R_init/ s/^/void rtmb_set_shared_pointers();\n/' RTMBconvenience/src/RcppExports.cpp
	sed -i '/R_useDynamicSymbols/ s/$$/\n    rtmb_set_shared_pointers();/' RTMBconvenience/src/RcppExports.cpp

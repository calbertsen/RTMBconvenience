R?=R-4-2-2

roxygen: 
	echo 'roxygen2::roxygenize("RTMBconvenience")' | $(R) --slave

clean:
	cd RTMBconvenience && ./cleanup

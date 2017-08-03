default:
	make clean
	mkdir -p results/current_run
	make debug
	make release
debug:
	bash run_example.sh debug

release:
	bash run_example.sh release

clean:
	rm -rf *.diag *.log *.aux *.out debug/ release/
	#rm -rf results/current_run
	rm -rf results/failed.tests results/results.summary
	rm -f  failed.tests results.summary
	rm -rf doc/figures/*-eps-converted-to.pdf
	rm -rf Time*
	rm -rf tmp*
	rm -rf Monodomain*

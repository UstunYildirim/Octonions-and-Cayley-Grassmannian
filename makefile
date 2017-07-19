main:
	ghc -O2 -outputdir aux/ Main.hs
mainThreaded:
	ghc -threaded -O2 -outputdir aux/ Main.hs
mainNoOptimization:
	ghc -outputdir aux/ Main.hs
profile:
	ghc --make -outputdir aux/ Main.hs -prof -auto-all -fforce-recomp
run-profile:
	./Main +RTS -p
clean:
	rm aux/*

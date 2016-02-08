all: Pos.txt.index

Pos.txt.index: Pos.txt
	rm -f Pos.txt.index
	rm -f Pos.txt.obd
	polyview Pos.txt

Pos.txt : main
	./main

PlotOrientierung.pdf: PlotLCOrient.py Orientierung.txt
	python PlotLCOrient.py

Orientierung.txt: main
	./main

VerteilungN_8.pdf: plotVerteilungen.py VerteilungN_8.txt
	python plotVerteilungen.py

VerteilungN_8.txt: main
	./main

PlotEndToEnd.pdf: PlotEndtoEnd.py EndToEnd.txt
	python PlotEndtoEnd.py

EndToEnd.txt: main
	./main

main : main.cpp bead.cpp bead.h vektor.cpp vektor.h chain.cpp chain.h global.cpp contactmap.h contactmap.cpp
	g++ -std=c++11 -O2 -lpthread main.cpp -o main

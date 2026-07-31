# compile cli tool
make tool

# find help regarding functions available in cli
./adapt -h
./adapt radar -h
./adapt viewshed -h
./adapt run -h

# generate viewshed (srinager airbase as source)
./adapt viewshed --lat 33.9856 --lon 74.7712 -a 200 -r 300 --tif
./adapt viewshed --lat 33.9856 --lon 74.7712 -a 200 -r 300 --svg
./adapt viewshed --lat 33.9856 --lon 74.7712 -a 200 -r 300 --png

./adapt viewshed --lat 33.9856 --lon 74.7712 -a 200 -r 300 --use-curvature --tif
./adapt viewshed --lat 33.9856 --lon 74.7712 -a 200 -r 300 --use-curvature --svg
./adapt viewshed --lat 33.9856 --lon 74.7712 -a 200 -r 300 --use-curvature --png

./adapt viewshed --lat 33.9856 --lon 74.7712 -a 200 -r 300 -s --loa 60 --roa 150 --tif
./adapt viewshed --lat 33.9856 --lon 74.7712 -a 200 -r 300 -s --loa 60 --roa 150 --svg
./adapt viewshed --lat 33.9856 --lon 74.7712 -a 200 -r 300 -s --loa 60 --roa 150 --png

# generate radar detection profile (srinager airbase as source)
./adapt radar --lat 33.9856 --lon 74.7712 -a 20 -t 2500 -r 300 --min-elev -5 --max-elev 60 --tif
./adapt radar --lat 33.9856 --lon 74.7712 -a 20 -t 2500 -r 300 --min-elev -5 --max-elev 60 --svg
./adapt radar --lat 33.9856 --lon 74.7712 -a 20 -t 2500 -r 300 --min-elev -5 --max-elev 60 --png

./adapt radar --lat 33.9856 --lon 74.7712 -a 20 -t 2500 -r 300 --min-elev -5 --max-elev 60 -s --loa 300 --roa 90 --tif
./adapt radar --lat 33.9856 --lon 74.7712 -a 20 -t 2500 -r 300 --min-elev -5 --max-elev 60 -s --loa 300 --roa 90 --svg
./adapt radar --lat 33.9856 --lon 74.7712 -a 20 -t 2500 -r 300 --min-elev -5 --max-elev 60 -s --loa 300 --roa 90 --png

./adapt radar --lat 33.9856 --lon 74.7712 -a 20 -t 2500 -r 300 --min-elev -5 --max-elev 60 -s --loa 300 --roa 90 --use-curvature --tif
./adapt radar --lat 33.9856 --lon 74.7712 -a 20 -t 2500 -r 300 --min-elev -5 --max-elev 60 -s --loa 300 --roa 90 --use-curvature --svg
./adapt radar --lat 33.9856 --lon 74.7712 -a 20 -t 2500 -r 300 --min-elev -5 --max-elev 60 -s --loa 300 --roa 90 --use-curvature --png

# generate radar detection in fast mode (srinager airbase as source)
./adapt radar --fast --lat 33.9856 --lon 74.7712 -a 20 -t 2500 -r 300 --min-elev -5 --max-elev 60 -s --loa 300 --roa 90 --tif

# generate near-optimal solutions using ga for homo, hetero and auto cases (update config accordingly)
./adapt run 
./adapt run --png

# benchmark ga strategies for homo, hetero and auto cases(update config accordingly)
./adapt run -b

# store solutions in h5 file for 100 trials (to develop ml model for heatmap generation)
./adapt run -g 5

# compile cli tool
make tool

# find help regarding functions available in cli
adapt -h
adapt radar -h
adapt viewshed -h
adapt run -h

# generate viewshed without sector (srinager airbase as source)
adapt viewshed --lat 33.9856 --lon 74.7712 -a 200 -r 300

# generate viewshed without sector (srinager airbase as source) - with curvature corrections
adapt viewshed --lat 33.9856 --lon 74.7712 -a 200 -r 300 --use-curvature

# generate viewshed with sector (srinager airbase as source)
adapt viewshed --lat 33.9856 --lon 74.7712 -a 200 -r 300 -s --loa 60 --roa 150

# generate radar detection profile without sector, with output image (srinager airbase as source)
adapt radar --lat 33.9856 --lon 74.7712 -a 20 -t 2500 -r 300 --min-elev -5 --max-elev 60

# generate radar detection profile without sector, without output image, in fast mode (srinager airbase as source)
adapt radar --fast --lat 33.9856 --lon 74.7712 -a 20 -t 2500 -r 300 --min-elev -5 --max-elev 60

# generate radar detection profile with sector, with output image  (srinager airbase as source)
adapt radar --lat 33.9856 --lon 74.7712 -a 20 -t 2500 -r 300 --min-elev -5 --max-elev 60 -s --loa 300 --roa 90

# generate radar detection profile with sector, without output image, in fast mode (srinager airbase as source)
adapt radar --fast --lat 33.9856 --lon 74.7712 -a 20 -t 2500 -r 300 --min-elev -5 --max-elev 60 -s --loa 300 --roa 90

# generate near-optimal solutions using ga for homo, hetero and auto cases (update config accordingly)
adapt run

# benchmark ga strategies for homo, hetero and auto cases(update config accordingly)
adapt run -b

# store solutions in h5 file for 100 trials (to develop ml model for heatmap generation)
adapt run -g 5

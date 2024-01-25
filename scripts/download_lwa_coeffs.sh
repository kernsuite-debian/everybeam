# Download LWA coefficients file
# These coefficients are generated using the script in 
# 'EveryBeam/scripts/coeff_scripts/convert_lwa.py', using the input simulation
# files available at https://git.astron.nl/RD/EveryBeam/-/issues/74#note_71706

set -e
LWA_COEFFICIENTS_FILE="LWA_OVRO.h5"

if [ ! -f "$LWA_COEFFICIENTS_FILE" ]; then
    wget https://support.astron.nl/software/ci_data/EveryBeam/${LWA_COEFFICIENTS_FILE}
fi
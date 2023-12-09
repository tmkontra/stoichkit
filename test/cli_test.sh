set -euxo pipefail

./target/debug/stoichkit balance H2O O2 = H2O2

./target/debug/stoichkit theoretical-yield "2*H2O2" 4.0 = "2*H2O" O2

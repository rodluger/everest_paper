cd figures

echo "Running cdpp.py..."
python cdpp.py

echo "Running injections.py..."
python injections.py

echo "Running nPLD.py..."
python nPLD.py

echo "Running outliers.py..."
python outliers.py

echo "Running ridge_reg.py..."
python ridge_reg.py

echo "Running saturated_star.py..."
python saturated_star.py
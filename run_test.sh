#/bin/bash

for ((i = 0; i < 100; i++)) do
    python3 discovery.py >> result.txt
done

python3 create_average.py
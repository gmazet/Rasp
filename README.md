# Rasberryshake
Few python code to extract and view Raspberryshake data

1. First set TMP environnement variable.
ex: export TMP=/tmp

2. Then execute: python ./syncRASPstations.py
This will download Raspberryshake station book

3. Edit mysta.py and put your own Rasp station codes
If you want to plot more than 3 stations, use option --nbsta


Usage:
Ex1: Retrieve Rasp data of the closest stations for a specific EMSC EventID (see https://www.emsc-csem.org/Earthquake/)
python main.py --evid 773032

Ex2: Retrieve Rasp data of yor own station for a specific EMSC EventID
python main.py --evid 773032 --mysta

Ex3: Add filtering
python main.py --evid 773032 --mysta --fmin 0.2 --fmax 5.0

Ex4: Show only one station
python main.py --evid 773032 --mysta --nbsta 1

Ex5: Show more signal (default is 120 sec)
python main.py --evid 773032 --mysta --len 180 --force

Ex6: Don't use EMSC EventID but a custom epicenter coordinates et origin time
python main.py --mysta --otime '2019-06-16T01:03:27.0Z' --lat 46.1 --lon 0.2 --mag 2.6 --depth 2

Ex7: Request RESIF data instead of Rasp
python main.py --provider resif --otime '2019-07-16T08:03:27.0Z' --lat 46.1 --lon 0.2 --mag 2.6 --depth 2



curl -o rasp.sta "https://fdsnws.raspberryshakedata.com/fdsnws/station/1/query?endafter=2019-09-26T00%3A00%3A00&network=AM&channel=SHZ&includerestricted=true&format=text&formatted=true&nodata=404"

wc -l rasp.sta

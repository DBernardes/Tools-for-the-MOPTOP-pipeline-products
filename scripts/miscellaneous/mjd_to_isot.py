from astropy.time import Time

time = "59863.782465"
t = Time(time, format="mjd")
print(t.to_value("isot"))

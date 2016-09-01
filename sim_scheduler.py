from datetime import datetime, timedelta


curdate = datetime(2017,12,31,0,0,0)
nt = 0
while nt < 365:
    d=random.randint(2,4)
    td = timedelta(d)
    curdate = curdate + td
    curdatestr = curdate.strftime("%Y-%m-%d")
    print "sim_night -d %s" % (curdatestr)

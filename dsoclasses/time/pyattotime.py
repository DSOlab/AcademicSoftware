import attotime
import datetime

def at2pt(at):
    """ Warning! This will cuse loss of precision.
        Translate an attotime instance to a native python datetime instance.
    """
    return datetime.datetime(at.year, at.month, at.day, at.hour, at.minute, at.second, at.microsecond)

def fsec2asec(fsec):
    isec = int(fsec) # integral seconds
    imsec = int((fsec - isec)*1e6)  # integral microseconds
    fnsec = float(fsec*1e9-imsec*1e3) # fractional nanoseconds
    assert abs(float(attotime.attotimedelta(seconds=isec, microseconds=imsec, nanoseconds=fnsec).total_nanoseconds()) - fsec*1e9) < 1e-1
    return attotime.attotimedelta(seconds=isec, microseconds=imsec, nanoseconds=fnsec)

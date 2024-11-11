import attotime
import datetime

def at2pt(at):
    """ Warning! This will cuse loss of precision.
        Translate an attotime instance to a native python datetime instance.
    """
    return datetime.datetime(at.year, at.month, at.day, at.hour, at.minute, at.second, at.microsecond)

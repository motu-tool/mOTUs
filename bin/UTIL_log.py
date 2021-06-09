import sys
import datetime

def get_timestamp():
    return '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())


def print_log(text):
    sys.stderr.write(get_timestamp()+" "+text+"\n")
    sys.stderr.flush()

def print_message(text):
    sys.stderr.write("   "+text+"\n")
    sys.stderr.flush()

def print_message_time(text):
    sys.stderr.write("  "+get_timestamp()+" "+text+"\n")
    sys.stderr.flush()

def print_error(text, exit = True):
    sys.stderr.write("   "+"Error: "+text+"\n")
    sys.stderr.flush()
    if exit:
        sys.exit(1)

def print_warning(text):
    sys.stderr.write("   "+"Warning: "+text+"\n")
    sys.stderr.flush()

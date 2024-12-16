import datetime
import random
import string


def random_time_string(name=None):
    current_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")  # formated time
    random_chars = ''.join(random.choices(string.ascii_letters + string.digits, k=10))  # random string
    # combine them
    if name:
        time_string = f"{current_time}__{name}__{random_chars}"
    else:
        time_string = f"{current_time}__{random_chars}"

    return time_string

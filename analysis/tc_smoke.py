import os
from tc_python import TCPython, LoggingPolicy

os.environ.setdefault('TC25B_HOME', '/opt/Thermo-Calc/2025b')

with TCPython(logging_policy=LoggingPolicy.NONE) as session:
    session.disable_caching()
    print('ok')

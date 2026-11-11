"""
CalibrateContTimer

Runs CalibrateCont every 0.1 s (fixed spacing between runs), like the MATLAB
timer('TimerFcn','CalibrateCont','Period',0.1,'ExecutionMode','fixedSpacing').

Run this file directly, or import it and call start()/stop():
    import CalibrateContTimer
    t = CalibrateContTimer.start()
    ...
    CalibrateContTimer.stop(t)
"""
import threading
import time

from CalibrateCont import CalibrateCont


class _FixedSpacingTimer(threading.Thread):
    def __init__(self, fcn, period):
        super().__init__(daemon=True)
        self.fcn = fcn
        self.period = period
        self._stop_event = threading.Event()

    def run(self):
        while not self._stop_event.is_set():
            try:
                self.fcn()
            except Exception as e:  # like a MATLAB timer, stop when the callback errors
                print(f'CalibrateCont error, timer stopped: {e}')
                break
            self._stop_event.wait(self.period)

    def stop(self):
        self._stop_event.set()


def start(period=0.1):
    Timer = _FixedSpacingTimer(CalibrateCont, period)
    Timer.start()
    return Timer


def stop(Timer):
    Timer.stop()


if __name__ == '__main__':
    Timer = start()
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        stop(Timer)

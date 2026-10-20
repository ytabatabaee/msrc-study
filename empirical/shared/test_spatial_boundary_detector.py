import numpy as np
from spatial_boundary_detector import segment_track, intervals_to_coordinates

def test_internal_boundary():
    y=np.r_[np.zeros(20),np.ones(12)*4,np.zeros(20)][:,None]
    s=segment_track(y,['x'],min_size=3,max_intervals=2)
    assert s.status=='internal' and s.intervals
def test_censored_and_no_change():
    y=np.r_[np.ones(10)*3,np.zeros(25)][:,None]
    assert segment_track(y,['x'],min_size=3,max_intervals=1).status in {'left','internal','no_detection'}
    assert segment_track(np.zeros((30,1)),['x'],min_size=3).status=='no_detection'
def test_missing_and_multiple():
    y=np.zeros((60,1));y[8:18]=3;y[35:45]=-2
    s=segment_track(y,['x'],min_size=3,max_intervals=2)
    assert len(s.intervals)>=1
    assert intervals_to_coordinates(s,[{'start':i,'end':i} for i in range(60)])

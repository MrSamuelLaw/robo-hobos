# code

def step_info(t,yout):

    # find max over steady state
    OS = ((max(yout)/yout[-1]) - 1) * 100

    # loop over and find the first time where it reaches 90% of steady state
    Tr = (t[next(i for i in range(0,len(yout)-1) if yout[i]>yout[-1]*.90)]-t[0])

    # loop over and find where the system reaches steady state
    if OS:
        # if no overshoot, calculate assuming approach from below
        Ts = (t[next(len(yout)-i for i in range(2,len(yout)-1) if abs(yout[-i]/yout[-1])>1.02)]-t[0])
    else:
        # if overshoot, calculate assuming approach from above
        Ts = (t[next(len(yout)-i for i in range(2,len(yout)-1) if abs(yout[-1]/yout[-i])>1.02)]-t[0])

    return OS, Tr, Ts, f"OS = {OS}\nTr = {Tr}\nTs = {Ts}"

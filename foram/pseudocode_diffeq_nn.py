import elvet

def equations(x,y,dy):
    y1,y2,y3,y4,y5,y6,y7 = y[0], y[1], y[2], y[3], y[4], y[5], y[6]
    dy1,dy2,dy3,dy4,dy5,dy6,dy7 = dy[0,0], dy[0,1], dy[0,2], dy[0,3],                                     
                                  dy[0,4], dy[0,5], dy[0,6]
    return [-1/x * dy1 - L1(y1,...y7), #assuming function L1 is known
            ....
            -1/x * dy7 - L7(y1,...y7)] #returns 7 equations
                        
# A and B are arrays of shape (7,)
# a and b are int with a,b > 0
# bcs is a list of size 14
bcs = [elvet.BC(b, lambda x, y, dy: y[0] - B[0]),
       ....
       elvet.BC(b, lambda x, y, dy: y[6] - B[6]),
       elvet.BC(a, lambda x, y, dy: dy[0,0] - A[0],
       ...        
       elvet.BC(a, lambda x, y, dy: dy[0,6] - A[6]
       ] 


domain = elvet.box((a, b, 100))

y = elvet.nn(1, 10, 7)

solver = elvet.solver(equations, bcs, domain, model=y, epochs=60000)
"""Linear response to the leading product-formula commutator, in units w=1."""
import numpy as np
from PIL import Image, ImageDraw
from pathlib import Path
import schwinger as model

def study(n, mass, steps):
    model.N, model.MASS = n, mass
    h, gates, vac, obs = model.model()
    # Restore left-to-right chronological order: later operators multiply on the left.
    previous = np.zeros_like(h)
    correction = np.zeros_like(h)
    for _, p, c in gates:
        term = c*p
        correction += -0.5j*(term @ previous - previous @ term)
        previous += term
    e, v = np.linalg.eigh(h)
    k = v.conj().T @ correction @ v
    initial = v.conj().T @ vac
    times = np.linspace(0, 2, 81)
    phases = np.exp(-1j*times[:,None]*e)
    exact = (phases*initial) @ v.T
    gap = e[:,None]-e[None,:]
    response = []
    for t, phase, state in zip(times, phases, exact):
        # Divided difference of exp(-itE), including degenerate eigenspaces.
        derivative = -1j*t*np.exp(-0.5j*t*(e[:,None]+e[None,:]))*np.sinc(t*gap/(2*np.pi))
        delta = v @ ((derivative*k) @ initial)
        response.append(2*t*obs @ (state.conj()*delta).real)
    response = np.array(response)
    psi = np.tile(vac,(len(times),1)).astype(complex)
    for _ in range(steps):
        for _, p, c in gates:
            indices = abs(p).argmax(axis=1)
            phase = p[np.arange(len(p)),indices]
            angle = times[:,None]*c/steps
            psi = np.cos(angle)*psi-1j*np.sin(angle)*phase*psi[:,indices]
    error = (abs(psi)**2-abs(exact)**2) @ obs.T
    # Reflection plus spin flip interchanges mirrored fluxes in the neutral sector.
    assert np.max(abs(response[:, 0])) < 1e-11
    assert np.max(abs(response[:, 1:] + response[:, :0:-1])) < 1e-11
    residual = error-response/steps
    return [n,mass,steps,abs(response[:,0]).max(),abs(response[:,1:]).max(),
            abs(error[:,0]).max(),abs(error[:,1:]).max(),abs(residual[:,1:]).max()]

if __name__ == "__main__":
    rows = [study(n,m,r) for n in (4,6) for m in (0,0.5,1) for r in (16,32,64)]
    model.write_csv("response.csv",["sites","mass","steps","density_linear_coefficient",
                    "flux_linear_coefficient","density_error","flux_error","flux_residual"],rows)
    print(np.array(rows))
    im = Image.new("RGB",(750,410),"white")
    d = ImageDraw.Draw(im)
    d.text((55,15),"Four sites, m/w = 0.5   |   Maximum absolute error",fill="black",font_size=21)
    data = [row for row in rows if row[:2]==[4,0.5]]
    def y(value):
        return 70-85*(np.log10(value)+1)
    for power in (1,2,3,4):
        yy=y(10**-power)
        d.line((65,yy,670,yy),fill="lightgray")
        d.text((10,yy-8),f"1e-{power}",fill="black",font_size=15)
    for idx,name,color in [(5,"density","#226161"),(6,"flux","#ad5521"),(7,"flux residual","#513c98")]:
        points=[(100+210*i,y(row[idx])) for i,row in enumerate(data)]
        d.line(points,fill=color,width=3)
        for x,yy in points:
            d.ellipse((x-4,yy-4,x+4,yy+4),fill=color)
        d.text((540,points[-1][1]-6),name,fill=color,font_size=16)
    for i,row in enumerate(data):
        d.text((90+210*i,354),str(row[2]),fill="black",font_size=16)
    d.text((240,380),"Product-formula steps",fill="black",font_size=17)
    im.save(Path(__file__).with_name("convergence.png"))


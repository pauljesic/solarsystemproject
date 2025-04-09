from astroquery.jplhorizons import Horizons
from vpython import sphere, vector, color, rate, canvas, label, mag, norm
from datetime import datetime, timedelta

AU = 1.496e11  # meters
DAY = 86400    # seconds
G = 6.67430e-11  # gravitational constant in m^3 kg^-1 s^-2
SCALE = 0.00000001 #global scaling factor for visualization

scene = canvas(title='To scale solar system')
scene.background = color.black




masses = {
    "Sun": 1.989e30,
    "Mercury": 3.301e23,
    "Venus": 4.867e24,
    "Earth": 5.972e24,
    "Mars": 6.417e23,
    "Jupiter": 1.898e27,
    "Saturn": 5.683e26,
    "Uranus": 8.681e25,
    "Neptune": 1.024e26
}

planet_ids = {
    "Mercury": "199",
    "Venus": "299",
    "Earth": "399",
    "Mars": "499",
    "Jupiter": "599",
    "Saturn": "699",
    "Uranus": "799",
    "Neptune": "899"
}

colors_dict = {
    "Sun" : color.yellow,
    "Mercury": color.gray(0.5),
    "Venus": color.orange,   
    "Earth": color.blue,
    "Mars": color.red,
    "Jupiter": color.orange,
    "Saturn": color.yellow,
    "Uranus": color.cyan,
    "Neptune": color.purple
}

radii = {
    "Sun" : 6.9634e7,
    "Mercury": 2.44e6,
    "Venus": 6.051e6,
    "Earth": 6.371e6,
    "Mars": 3.389e6,
    "Jupiter": 6.9911e7,
    "Saturn": 5.8232e7,
    "Uranus": 2.5362e7,
    "Neptune": 2.4622e7
}

class Planet:
    def __init__(self, name, mass, position, velocity, radius, color):
        self.name = name
        self.mass = mass    #kg
        self.pos = position  #meters
        self.vel = velocity  #m/s
        self.radius = radius  #meters
        self.color = color

        self.sphere = sphere(
            pos=self.pos * SCALE,
            radius=max(self.radius * 0.0000015, 0.2),
            color=self.color,
            make_trail=True
        )
        self.label = label(pos=self.sphere.pos + vector(self.sphere.radius * 2, 0, 0), text=self.name, height=10)

    def update_visual(self):
        self.sphere.pos = self.pos * SCALE     #self.pos is not scaled so you have to scale it before passing 
        self.label.pos = self.sphere.pos + vector(self.sphere.radius * 2, 0, 0)



def get_planet_data(planet_id: str, date: str) -> dict:
    """
    retrieves position and velocity data from NASA JPL Horizons for a given planet and date

    planet_id (str): horizons ID for the planet ("399" for earth)
    date (str): The date in "YYYY-MM-DD" 

    Returns a dictionary with 'position' and 'velocity' lists
    """

    start_time = date
    #stop time is required because jpl requires a date a range
    stop_time = (datetime.strptime(date, "%Y-%m-%d %H:%M") + timedelta(days=1)).strftime("%Y-%m-%d %H:%M")

    #"step": "1d" ensure only 1 data point is taken from the date range. this is done because it is a lot simpler than setting up the julian time system that is otherwise requried
    obj = Horizons(id=planet_id, location='@sun', epochs={"start": start_time, "stop": stop_time, "step": "1d"})
    vectors = obj.vectors()

    return {
        'position': [
            float(vectors['x'][0]),
            float(vectors['y'][0]),
            float(vectors['z'][0])
        ],
        'velocity': [
            float(vectors['vx'][0]),
            float(vectors['vy'][0]),
            float(vectors['vz'][0])
        ]
    }




def compute_gravitational_force(body1, body2):
    r_vec = body2.pos - body1.pos       #get 3d cordinates between the 2 bodies
    distance = mag(r_vec)       #get actual distance between the bodies: sqrt(x^2 + y^2 + z^2)
    
    force_magnitude = G * body1.mass * body2.mass / distance**2     #GMm/r^2 gravititional force in N
    force_direction = norm(r_vec)       #gets direction of force 
    return force_magnitude * force_direction




if __name__ == "__main__":
    #astroquery requires hours and mintues so 00:00 is appended to the user inputted date
    date_input = input("Enter a date (YYYY-MM-DD): ") + " 00:00"
    

    #add the Sun as a Planet object with fixed position and zero velocity
    sun_obj = Planet(
        name="Sun",
        mass=masses["Sun"],
        position=vector(0, 0, 0),
        velocity=vector(0, 0, 0),
        radius=radii["Sun"],
        color=colors_dict["Sun"]
    )

    planet_objects = [sun_obj]

    for name, pid in planet_ids.items():
        print(f"Fetching data for {name}...")
        data = get_planet_data(pid, date_input)
        pos_m = vector(*data['position']) * AU # unpacks the list and turns it into a VPython 3D vector. horizons gives the value as AU so u have to multiply by AU
        vel_mps = vector(*data['velocity']) * AU / DAY  # converts AU/day (jpl horizons format) to m/s (so it can work with gravity formula)

        planet = Planet(
            name=name,
            mass=masses[name],
            position=pos_m,
            velocity=vel_mps,
            radius=radii[name],
            color=colors_dict[name]
        )
        planet_objects.append(planet)

    

    while True:
        rate(60)
        TIME_SCALE = 24  #simulate 24 hours per real second
        dt = 60 * 60 * TIME_SCALE  #one day per second

        #initiliaze empty dictionary with each key being a planet in planet_objects and vector velocity =0 
        forces = {planet: vector(0, 0, 0) for planet in planet_objects}

        #compute net force on each planet

        for i, planet in enumerate(planet_objects):     #nested loop: for every planet in planetobjs, go through each planet and calculate the forces acting on it
            for j, other in enumerate(planet_objects):  
                if i != j:                              #dont calculate the forces a planet is acting on itself
                    force = compute_gravitational_force(planet, other)      #calculate the planet acts on another planet
                    forces[planet] += force                 #adds total forces that each planet is acting on each other to a dictionary

        #update velocity and position using Euler integration
        for planet in planet_objects:
            if planet.name != "Sun":  # Keep the Sun stationary
                acceleration = forces[planet] / planet.mass #acell = force/mass
                planet.vel += acceleration * dt              #Updates the velocity of the planet based on the acceleration and time step dt  eulers intergration Vnew = Vold + a * dt
                planet.pos += planet.vel * dt                #Updates the position of the planet using its updated velocity     Rnew = Rold + v * dt
            planet.update_visual()

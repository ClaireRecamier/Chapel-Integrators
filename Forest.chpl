//define constants pi, solarMass, number of timesteps

config const iterations = 100000, timestep = 0.5;
const pi = 3.141592653589793, daysPerYear = 365.24, solarMass = 4 * pi**2;

//define class/record of a body
record body {
  var position: 3*real;
  var velocity: 3*real;
  var mass: real;
}

//define an array of bodies
var bodies = [
  //SUN
  new body(mass = solarMass),

  /* jupiter */
  new body(position = ( 4.84143144246472090e+00,
                  -1.16032004402742839e+00,
                  -1.03622044471123109e-01),
           velocity = ( 1.66007664274403694e-03 * daysPerYear,
                   7.69901118419740425e-03 * daysPerYear,
                  -6.90460016972063023e-05 * daysPerYear),
          mass =   9.54791938424326609e-04 * solarMass),

  /* saturn */
  new body(position = ( 8.34336671824457987e+00,
                   4.12479856412430479e+00,
                  -4.03523417114321381e-01),
           velocity = (-2.76742510726862411e-03 * daysPerYear,
                   4.99852801234917238e-03 * daysPerYear,
                   2.30417297573763929e-05 * daysPerYear),
          mass =   2.85885980666130812e-04 * solarMass),

  /* uranus */
  new body(position = ( 1.28943695621391310e+01,
                  -1.51111514016986312e+01,
                  -2.23307578892655734e-01),
           velocity = ( 2.96460137564761618e-03 * daysPerYear,
                   2.37847173959480950e-03 * daysPerYear,
                  -2.96589568540237556e-05 * daysPerYear),
          mass =   4.36624404335156298e-05 * solarMass),

  /* neptune */
  new body(position = ( 1.53796971148509165e+01,
                  -2.59193146099879641e+01,
                   1.79258772950371181e-01),
           velocity = ( 2.68067772490389322e-03 * daysPerYear,
                   1.62824170038242295e-03 * daysPerYear,
                  -9.51592254519715870e-05 * daysPerYear),
          mass =   5.15138902046611451e-05 * solarMass)
];
const numBodies = bodies.size;

//procedure to initialize velocity of sun. Checked, function works
proc initSun() {
  var p: 3*real = (0.0,0.0,0.0);
  for body in bodies do { //shortcut: could use reduce instead to shorten this function
    p += body.mass * body.velocity;
  }
  //const p = + reduce (for b in bodies do (b.velocity * b.mass));
  bodies[0].velocity = -p / solarMass;
  writeln("vel of sun:", bodies[0].velocity);
}

//procedure to advance by a timestep
proc advance(dt: real){
  //Forest integrator constants
  const a = 1.35120;
  const b = -0.3512;
  const c = -1.7024;

  //define force vector
  var F: 3*real = (0.0,0.0,0.0);
  //define array of Force tuples;
  var Forces: [0..numBodies-1] 3*real;
  //initialize relative position vector and its length as new variables
  var r: 3*real;
  var r_len: real = 0.0;

  //Drift by a*dt
  for i in 0..numBodies-1 do{
    bodies[i].position += a * dt * bodies[i].velocity;
  }

  //Kick by a*dt, drift by b*dt
  for i in 0..numBodies-1 do{
    for j in i+1..numBodies-1 do {
      //calculate relative position vector
      r = bodies[i].position - bodies[j].position;
      //calculate length of relative position vector
      r_len = dist(r);
      //Force btwn two bodies: F = G m1 m2 / r ^2 * unit vector
      F = ((bodies[i].mass * bodies[j].mass)/ r_len**2 ) * (r/r_len);
      Forces[i] -= F; //add this force to total force on body i
      Forces[j] += F; //add this force to total force on body j
    }
    //add acceleration vector to velocity vector for new velocity
    bodies[i].velocity += dt * Forces[i]/bodies[i].mass * a;
    bodies[i].position += dt * bodies[i].velocity * b;
  }

  //Kick by c*dt, drift by b*dt
  for i in 0..numBodies-1 do{
    for j in i+1..numBodies-1 do {
      //calculate relative position vector
      r = bodies[i].position - bodies[j].position;
      //calculate length of relative position vector
      r_len = dist(r);
      //Force btwn two bodies: F = G m1 m2 / r ^2 * unit vector
      F = ((bodies[i].mass * bodies[j].mass)/ r_len**2 ) * (r/r_len);
      Forces[i] -= F; //add this force to total force on body i
      Forces[j] += F; //add this force to total force on body j
    }
    //add acceleration vector to velocity vector for new velocity
    bodies[i].velocity += dt * Forces[i]/bodies[i].mass * c;
    bodies[i].position += dt * bodies[i].velocity * b;
  }
  /*
  //Kick by a*dt, drift by a*dt
  for i in 0..numBodies-1 do{
    for j in i+1..numBodies-1 do {
      //calculate relative position vector
      r = bodies[i].position - bodies[j].position;
      //calculate length of relative position vector
      r_len = dist(r);
      //Force btwn two bodies: F = G m1 m2 / r ^2 * unit vector
      F = ((bodies[i].mass * bodies[j].mass)/ r_len**2 ) * (r/r_len);
      Forces[i] -= F; //add this force to total force on body i
      Forces[j] += F; //add this force to total force on body j
    }
    //add acceleration vector to velocity vector for new velocity
    bodies[i].velocity += dt * Forces[i]/bodies[i].mass * a;
    bodies[i].position += dt * bodies[i].velocity * a;
  }
  */


}

//procedure to calculate energy
proc energy(){
  var energy: real = 0.0;
  //initialize relative position vector and its length as new variables
  var r: 3*real = (0.0,0.0,0.0);
  var r_len: real = 0.0;

  for i in 0..numBodies-1 do {
    //add KE
    energy += 0.5 * bodies[i].mass * dist(bodies[i].velocity)**2;
    //subtract PE of interaction of this body with every other body (of which PE between has not already been calculated)
    for j in i+1..numBodies-1 do {
      //calculate relative position vector
      r = bodies[i].position - bodies[j].position;
      //calculate length of relative position vector
      r_len = dist(r);

      energy -= bodies[i].mass * bodies[j].mass / r_len;
    }
  }
  return energy;


}

//shortcut procedure to help calculate distance. Checked, function works.
proc dist(r) { //shortcut: pass in array as (x,y,z) instead.
  var x: real = r(0);
  var y: real = r(1);
  var z: real = r(2);
  var dist: real = sqrt(x**2 + y**2 + z**2);
  return dist;
}


//write main procedure. Initialize velocity of sun, print energy. Advance by a timestep, print energy.
proc main() {

  initSun();
  writef("%.9r\n",energy());

  for i in 1..iterations do {
    advance(timestep);
    if i % 1000 == 0 { //print energy every 100th iteration
      writef("%.9r\n",energy());
    }
  }
  //writef("%.9r\n",energy());


  /* alternative main proc, where change in energy prints instead of energy
  initSun();

  var init_energy: real = energy();
  var curr_energy: real = init_energy;
  var delta_energy: real = 0.0;
  writef("%.9r\n",init_energy);

  for i in 1..iterations do {
    advance(timestep);
    if i % 1000 == 0 { //print energy every 100th iteration
      curr_energy = energy();
      delta_energy = (curr_energy - init_energy)/init_energy;
      writef("%.9r\n",abs(delta_energy));
    }
  }
  //writef("%.9r\n",energy());
  */

}

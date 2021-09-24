#version 330 core

#define PI 3.14159265358979

#define BACKGROUND_COLOR vec3(0.5,0.8,0.9)//vec3(0.6,0.8,1.0)
#define FOG_COLOR vec3(0.30, 0.36, 0.40)
#define LIGHT_COLOR vec3(1.0,0.9,0.7)
#define LIGHT_DIRECTION vec3(0.36, 0.48, 0.80)

#define MAX_DIST 200.
#define MIN_DIST 1e-5
#define MAX_MARCHES 1000
#define SUN_SIZE 0.001
#define SUN_SHARPNESS 1.5
#define POWER 8
#define OSCILLATION 1
#define AA 1


in vec4 gl_FragCoord;
in vec2 gl_PointCoord;

out vec4 FragColor; 
uniform vec3 cameraPos;
uniform float iTime;
uniform vec2 iResolution;
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;


float res;

// https://github.com/HackerPoet/MarbleMarcher/blob/master/assets/frag.glsl
vec3 boxFold(vec3 z, vec3 r) {
	return clamp(z.xyz, -r, r) * 2.0 - z.xyz;
}

// http://www.fractalforums.com/fragmentarium/fragmentarium-an-ide-for-exploring-3d-fractals-and-other-systems-on-the-gpu/15/
void sphereFold(inout vec3 z, inout float dz) {
    float fixedRadius2 = .6 + 4.* cos(20./8.) + 4.;
    float minRadius2 = 0.3;
	float r2 = dot(z,z);
	if (r2< minRadius2) {
		float temp = (fixedRadius2/minRadius2);
		z*= temp;
		dz*=temp;
	} 
    else if (r2<fixedRadius2) {
		float temp =(fixedRadius2/r2);
		z*=temp;
		dz*=temp;
	}
}

// https://github.com/HackerPoet/MarbleMarcher/blob/master/assets/frag.glsl
vec3 mengerFold(vec3 z) {
	float a = max(z.x - z.y, 0.0);//min
	z.x -= a;
	z.y += a;
	a = max(z.x - z.z, 0.0);//min
	z.x -= a;
	z.z += a;
	a = max(z.y - z.z, 0.0);//min
	z.y -= a;
	z.z += a;
    return z;
}

float mandelbox(vec3 z)
{
    float scale = 2.0;
	vec3 offset = z;
	float dr = 1.0;
	for (int n = 0; n < 10; n++)
    {
		z = boxFold(z,vec3(2.0)); //CHANGE VEC3 TO CHANGE THE FRACTAL
		sphereFold(z,dr);
        z = scale * z + offset;
        dr = dr * abs(scale) + 1.0;
	}
	float r = length(z);
	return r / abs(dr);
}

vec2 mandelBox2(vec3 z)
{
    float Iterations = 20.;
    float Scale = 3;
	vec3 offset = z;
	float dr = 1.0;
    float trap = 1e10;
	for (float n = 0.; n < Iterations; n++) {      
        z = mengerFold(z);
        z = boxFold(z, vec3(2.));       // Reflect
        sphereFold(z, dr);    // Sphere Inversion
        z.xz = -z.zx;
		z = boxFold(z, vec3(1.));       // Reflect
        
		sphereFold(z, dr);    // Sphere Inversion
        z=Scale*z + offset;  // Scale & Translate
        dr = dr*abs(Scale)+1.0;
        trap = min(trap, length(z));
	}
	float r = length(z);
	return vec2(r/abs(dr), trap);
}

float sdMandelbulb(vec3 p) {
	vec3 z = p;
    float power = POWER;
    if(OSCILLATION == 1)
        power = (sin(iTime*0.1)+1) * POWER/2 + 1;
	float dr = 1.0;
	float r = 0.0;
	for (int i = 0; i < 15 ; i++) {
		r = length(z);
		if (r>2) break;	
		// convert to polar coordinates
		float theta = power*acos(z.z/r);
		float phi = power*atan(z.y,z.x);
		dr =  power*pow( r, power-1)*dr + 1.0;	
		// scale and rotate the point	
		// convert back to cartesian coordinates
		z = pow(r,power)*vec3(sin(theta)*cos(phi), sin(phi)*sin(theta), cos(theta));
		z+=p;
	}
	return 0.5*log(r)*r/dr;
}

float sdMandelbulb2(vec3 pos) {
	vec3 z = pos;
	float dr = 1.0;
	float r = 0.0;
	for (int i = 0; i < 15 ; i++) {
		r = length(z);
		if (r>2) break;	
		// convert to polar coordinates
		float theta = asin( z.z/r );
        float phi = atan( z.y,z.x );
		dr =  pow( r, POWER-1)*POWER*dr + 1.0;	
		// scale and rotate the point
		float zr = pow( r,POWER);
		theta = theta*POWER;
		phi = phi*POWER;		
		// convert back to cartesian coordinates
		z = zr*vec3( cos(theta)*cos(phi), cos(theta)*sin(phi), sin(theta) );
		z+=pos;
	}
	return 0.5*log(r)*r/dr;
}

float sdSphere(vec3 p, float s){
    return length(p)-s;
}

float sdTorus( vec3 p, vec2 t ){
  vec2 q = vec2(length(p.xz)-t.x,p.y);
  return length(q)-t.y;
}

float sdSphereMod(vec3 p, float s){
    vec3 sphere = vec3 (1.0,1.0,1.0);
    return length(mod(sphere.xyz -p,s) - vec3(s/2.0)) - .5;
    return length(p)-s;
}

float sdBox(vec3 p, vec3 b){
  vec3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

float sdMenger(vec3 p){
    float size=10.;
	p.z -=3.; 
    vec3[] s = vec3[](vec3(1,1,1),vec3(1,1,0));
    
    for(int iter=0;iter<10;++iter){
       
        p=abs(p);
        if(p.y > p.x) p.yx = p.xy;
        if(p.z > p.y) p.zy = p.yz;
        
        if(p.z > .5*size) p -= size*s[0];
        else p -= size*s[1];
        size /=3.;
        
    }
    float result = sdBox(p,vec3(1.5*size));
    return result;
}

float sdPlane(vec3 p){
    return p.y;
} 

float sminCubic( float a, float b, float k ){
    float h = max( k-abs(a-b), 0.0 )/k;
    return min( a, b ) - h*h*h*k*(1.0/6.0);
}

float distanceField(vec3 p){
    
    //float Sphere = sdSphere(p-vec3(0.0,0.0,0.0),2.0);
    //float Box = sdBox(p,vec3(2.0));
    //float Torus = sdTorus(p,vec2(2.0,0.5));
    //float SphereMod = sdSphereMod(p,2.0);
    //float Plane = sdPlane(p-vec3(0.0,-2.0,0.0));
    
    //float Menger = sdMenger(p);
    //float Mandelbulb2 = sdMandelbulb2(p);
    //float MandelBox = mandelbox(p);
    //float MandelBox2 = mandelBox2(p).x;
    //return sminCubic(Sphere, Torus,0.5);
    //return max(Torus,SphereMod);
    //return min(Sphere,Plane);
    
    float Box = sdBox(p,vec3(1.5));
    if(Box > 0.1) return Box;
    float Mandelbulb = sdMandelbulb(p);
    return Mandelbulb;

    //return SphereMod;
}

vec3 calcNormal(vec3 p, float h){ // https://www.iquilezles.org/www/articles/normalsSDF/normalsSDF.htm
    vec3 k = vec3(1,-1,0);
    return normalize( k.xyy*distanceField( p + k.xyy *h) + 
                      k.yyx*distanceField( p + k.yyx *h) + 
                      k.yxy*distanceField( p + k.yxy *h) + 
                      k.xxx*distanceField( p + k.xxx *h) );
}

//TO WORK THE BEST THIS NEED A PROPER DISTANCE FIELD COULD CAUSE ARTIFACTS WITH MIN()
float softShadow(in vec3 ro, in vec3 rd, float mint, float maxt, float k, float min_dist){ //https://www.iquilezles.org/www/articles/rmshadows/rmshadows.htm
    float res = 1.0;
    float ph = 1e10;
    for(float t=mint; t<maxt;){
        float h = distanceField(ro + rd*t);
        if( h < min_dist)
            return 0.0;
        float y = h*h/(2.0*ph);
        float d = sqrt(h*h-y*y);
        res = min(res, k*d/max(0.0,t-y));
        ph = h;
        t += h;  
    }
    res = clamp( res, 0.0, 1.0 );
    return res*res*(3.0-2.0*res);
    //return res;
}

float ambientOcclusion(vec3 p, vec3 n){
    float steps = 0.1;
    float ao = 0.0;
    float dist;
    for(int i=1; i<= 3;i++){
        dist = steps * i;
        ao += max(0.0,(dist - distanceField(p+n*dist))/ dist);
    }
    return (1.0-ao * 0.22);
}


vec3 raymarching(vec3 ro, vec3 rd, out int iter){
    float min_d = 1.0;
    float t = 0; //distance traveled alongside the ray vector
    float s = 0;
    float d = 0;
    for(int i = 0; i < MAX_MARCHES; i++){   
        iter = i;    
        float min_dist = max(res*t, MIN_DIST);
        //float min_d = min(min_d, 10.0 * d / t); Check how close we got without hitting so we can use it for something eg GLOW EFFECT
        d = distanceField(ro + rd * t);

        if (t > MAX_DIST){ break; }
        else if (d < min_dist){
            s += d / min_dist;         // can use this for fake ambient occlusion
            break;
        }     
        t += d;      
    }
    return vec3(d,t,s);
} 

vec4 render(vec3 ro, vec3 rd){
    vec4 col = vec4(0.0);
    int iter;
    vec3 raymarch = raymarching(ro,rd,iter);
    float d = raymarch.x;
    float t = raymarch.y;
    float s = raymarch.z;

    float min_dist = max(res*t, MIN_DIST);
    vec3 p = ro + rd * t;
    if(d < min_dist){
        vec3 n = calcNormal(p,min_dist);

        vec3 sun_light = LIGHT_COLOR * clamp(dot(n,LIGHT_DIRECTION ), 0., 1.);//LIGHT_COLOR * max(dot(n, LIGHT_DIRECTION), 0.0);
        vec3 sky_light = (BACKGROUND_COLOR*0.10)* clamp(0.5+0.5*dot(n,vec3(0.,1.,0.)), 0., 1.);
        vec3 bounce_light = (vec3(.06,.063,.07))* clamp(0.5+0.5*dot(n,vec3(0.,-1.,0.)), 0., 1.);
        //vec3 bounce_light = (LIGHT_COLOR*0.015)* sqrt(clamp( 0.1-0.9*n.y, 0.0, 1.0 ))*clamp(1.0-0.1*p.y,0.0,1.0);
        vec3 sun_shadow = LIGHT_COLOR * softShadow(p+n*min_dist, LIGHT_DIRECTION, 0.01, 64.0, 12.0, min_dist);
        
        // IF YOU HAVE A PLANE IT SHOULD TECNICALLY BOUNCE SOME LIGHT FROM THE SUN UNDER YOUR OBJECT, THIS LIGHT SHOULD
        // COME FROM THE BOTTOM dot(n,vec3(0.0,-1.0,0.0)) AND IT SHOULD BE THE COLOR OF THE SUN * THE PLANE OR BOTTOM LIGHT

        float ao = ambientOcclusion(p,n);
        
        /*float k = min_d * min(t, 1.0);
        vec3 reflected = n - 2.0 * dot(n,n) * n;
        float specular = max(dot(reflected, LIGHT_DIRECTION), 0.0);
        specular = pow(specular, 40);   
        
        //col.xyz += LIGHT_COLOR * specular * (0.25); //LIGHT COLOR * SPECULAR INTENSITY * SPECULAR FACTOR*/
        col.xyz += sun_light * sun_shadow;        
        col.xyz += sky_light * ao;   
        col.xyz += bounce_light * ao;
    }
    else {
        vec3 sky = BACKGROUND_COLOR - max(rd.y,0.0)*0.5;
        float sun = dot(rd, LIGHT_DIRECTION) - 1.0 + SUN_SIZE;
        sun = min(exp(sun * SUN_SHARPNESS / SUN_SIZE), 1.0);
        col.xyz += sky;
        col.xyz += LIGHT_COLOR * sun;
        //float fogAmount = 1.0 - exp(-(p.z-8.0) * (1.0/MAX_DIST));
        //col.xyz += mix(col.xyz, FOG_COLOR, fogAmount);
    }
    //return vec4(vec3(float(iter)/128.0), 1.0); //check ray iteration
    return col;
} 

void main()
{   
    vec4 tot;
    res = 1.0 / 2160.0;
    vec4 col = vec4(0.0);
    for (int i = 0; i < AA; ++i) {
		for (int j = 0; j < AA; ++j) {
            vec2 aa =-0.5+vec2(i,j)/AA;
            //vec2 uv = (gl_FragCoord.xy+0.5*(-iResolution.xy))/iResolution.y;
            vec2 uv = (((gl_FragCoord.xy+aa)-0.5*iResolution.xy))/iResolution.y; //vec2 uv = (-iResolution.xy + (2*gl_FragCoord.xy+aa))/iResolution.y;
            //AA cause object to distort on the edge of the screen to compensate you increase the FOCAL DISTANCE
           
            vec4 far_4 = inverse(view) * vec4(normalize(vec3(uv.x,uv.y,-1.41)), 1.0);

            vec3 far_3 = far_4.xyz/far_4.w;

            vec3 rd = normalize(far_3.xyz - cameraPos);

            //col.xyz =  vec3(length(uv),length(uv),length(uv));
            col = render(cameraPos,rd);

            col = pow(col, vec4(0.4545)); //GAMMA CORRECTION
            tot += col;                   
        }
    }
    tot /= AA*AA;
    FragColor = tot;
}
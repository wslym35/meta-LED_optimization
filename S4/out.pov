// -w320 -h240

#version 3.6;

#include "colors.inc"
#include "textures.inc"
#include "shapes.inc"

global_settings {max_trace_level 5 assumed_gamma 1.0}

camera {
	location <-1.620000, 3.240000, -4.860000>
	direction <0, 0,  2.25>
	right x*1.33
	look_at <0,0,0>
}

#declare Dist=80.0;
light_source {< -25, 50, -50> color White
	fade_distance Dist fade_power 2
}
light_source {< 50, 10,  -4> color Gray30
	fade_distance Dist fade_power 2
}
light_source {< 0, 100,  0> color Gray30
	fade_distance Dist fade_power 2
}

sky_sphere {
	pigment {
		gradient y
		color_map {
			[0, 1  color White color White]
		}
	}
}

#declare Xaxis = union{
	cylinder{
		<0,0,0>,<0.8,0,0>,0.05
	}
	cone{
		<0.8,0,0>, 0.1, <1,0,0>, 0
	}
	texture { pigment { color Red } }
}
#declare Yaxis = union{
	cylinder{
		<0,0,0>,<0,0.8,0>,0.05
	}
	cone{
		<0,0.8,0>, 0.1, <0,1,0>, 0
	}
	texture { pigment { color Green } }
}
#declare Zaxis = union{
	cylinder{
	<0,0,0>,<0,0,0.8>,0.05
	}
	cone{
		<0,0,0.8>, 0.1, <0,0,1>, 0
	}
	texture { pigment { color Blue } }
}
#declare Axes = union{
	object { Xaxis }
	object { Yaxis }
	object { Zaxis }
}
#declare Material_vac = texture{ pigment{ rgb <0.952230,0.513401,0.364784> } }
#declare Material_c-GaN(475nm) = texture{ pigment{ rgb <0.717297,0.635712,0.916195> } }
#declare Material_ITO(475nm) = texture{ pigment{ rgb <0.016301,0.606969,0.141603> } }
#declare Material_c-Sapp(475nm) = texture{ pigment{ rgb <0.804177,0.137232,0.242887> } }
#declare Material_c-GaN(540nm) = texture{ pigment{ rgb <0.129790,0.400944,0.156679> } }
#declare Material_c-Sapp(540nm) = texture{ pigment{ rgb <0.218257,0.998925,0.108809> } }
#declare Layer_sapp = union{
	difference{
		intersection{
			plane{ <0.540000,0.000000,0>, 0.270000 }
			plane{ <-0.540000,-0.000000,0>, 0.270000 }
			plane{ <0.000000,0.540000,0>, 0.270000 }
			plane{ <-0.000000,-0.540000,0>, 0.270000 }
			plane{ <0.540000,0.540000,0>, 0.381838 }
			plane{ <-0.540000,-0.540000,0>, 0.381838 }
			plane{ <0,0,-1>, 0 }
			plane{ <0,0,1>, 0.000000 }
		}
// nshapes = 0
		texture { Material_c-sapp(540nm) }
	}
	translate +z*0.000000
}
#declare Layer_GaN = union{
	difference{
		intersection{
			plane{ <0.540000,0.000000,0>, 0.270000 }
			plane{ <-0.540000,-0.000000,0>, 0.270000 }
			plane{ <0.000000,0.540000,0>, 0.270000 }
			plane{ <-0.000000,-0.540000,0>, 0.270000 }
			plane{ <0.540000,0.540000,0>, 0.381838 }
			plane{ <-0.540000,-0.540000,0>, 0.381838 }
			plane{ <0,0,-1>, 0 }
			plane{ <0,0,1>, 0.996000 }
		}
// nshapes = 2
box{
	<-1,-1,0>, <1,1,0.996000>
	scale +x*0.190000
	scale +y*0.270000
	rotate +z*0.000000
	translate +x*0.410000
	translate +y*0.000000
}
box{
	<-1,-1,0>, <1,1,0.996000>
	scale +x*0.145000
	scale +y*0.270000
	rotate +z*0.000000
	translate +x*0.133000
	translate +y*0.000000
}
		texture { Material_vac }
	}
	difference{
		intersection{
box{
	<-1,-1,0>, <1,1,0.996000>
	scale +x*0.190000
	scale +y*0.270000
	rotate +z*0.000000
	translate +x*0.410000
	translate +y*0.000000
}
			plane{ <0,0,-1>, 0 }
			plane{ <0,0,1>, 0.996000 }
		}
		texture { Material_c-GaN(540nm) }
	}
	difference{
		intersection{
box{
	<-1,-1,0>, <1,1,0.996000>
	scale +x*0.145000
	scale +y*0.270000
	rotate +z*0.000000
	translate +x*0.133000
	translate +y*0.000000
}
			plane{ <0,0,-1>, 0 }
			plane{ <0,0,1>, 0.996000 }
		}
		texture { Material_c-GaN(540nm) }
	}
	translate +z*0.000000
}
#declare Layer_air = union{
	difference{
		intersection{
			plane{ <0.540000,0.000000,0>, 0.270000 }
			plane{ <-0.540000,-0.000000,0>, 0.270000 }
			plane{ <0.000000,0.540000,0>, 0.270000 }
			plane{ <-0.000000,-0.540000,0>, 0.270000 }
			plane{ <0.540000,0.540000,0>, 0.381838 }
			plane{ <-0.540000,-0.540000,0>, 0.381838 }
			plane{ <0,0,-1>, 0 }
			plane{ <0,0,1>, 0.000000 }
		}
// nshapes = 0
		texture { Material_vac }
	}
	translate +z*0.996000
}
#declare Layers = union {
	//object{ Layer_sapp }
	object{ Layer_GaN }
	//object{ Layer_air }
}

Axes
Layers

# smallpt

Explicit lighting version of smallpt

## Environment:

Linux Mint

## Commands:

* ```make; ./pt_cpu glut 4 null```

* ./pt_cpu {"glut"|"console"} {samplePerPixels} {objPath}

* pre-setting scenes: 
    
    * null: cornell box

    * bunny: stanford bunnies

    * haku: yowane haku

## Features:

* Supports Polygon Models (.obj)

* Supports Polygon Area Lights

* OpenMp (CPU) Accelerated

* GLUT based window display current results

## TODO:

* Sample texture texels is in progress, comments are leave in codes

## Results:

* Original Scenes

![](Demos/125589.png)

* Soft Shadows

![](Demos/s_output_bunnyLowWithUV_63spp.png)

![](Demos/b_output_bunnyLowWithUV_63spp.png)

* Hard Shadows

![](Demos/h_output_bunnyLowWithUV_63spp.png)

* Polygon Area Lights

![](Demos/output_bunnyLowWithUV_32spp.png)

![](Demos/output_bunnyLowWithUV_16spp.png)

* Misc

![](Demos/output_b_5spp.png)

![](Demos/output_bunnyLowWithUV_17spp.png)

![](Demos/haku_output_2spp.png)
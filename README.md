# Visibility Polygon

Simple 2D Field of View/Line of Sight algorithm implemented in gdscript (the Godot Engine Scripting language). 

![visualization of simple fov examples](https://github.com/fahall/godot_2d_visibility/blob/master/fov_example.png "Simple Examples")

## Algorithm
This implementation is based on [this 2d Visibility](http://www.redblobgames.com/articles/visibility/) article on the Red Blob Games website. It also draws heavily from [Trylock's C++ implementation](https://github.com/trylock/visibility)


## Preconditions:
- The line segments must not intersect except at their endpoints 
- The visiblity polygon has to be closed. 

**Behaviour is undefined if the preconditions aren't met**. Line segment intersections can be solved by first splitting segments at intersection points. The visibility polygon can be guaranteed closed if you add an additional bounding box around your scene. 

## Sample Use in Godot:

This sample will draw a blue polygon indicating the computed FOV from the node's origin. 

All occluding objects are assumed to have a `polygon node` called `occlusion_mask` and be a member of the `fov_occluders` group. 

```
func _draw():
	fov()

func get_occluder_polygons():
	var occluders = get_tree().get_nodes_in_group('fov_occluders')
	var polygons = []
	for o in occluders:
		var mask = o.get_node('occlusion_mask')
		if mask != null:
			var poly = mask.polygon
			for i in range(len(poly)):
				poly[i] += (o.position - self.position)
			polygons.append(poly)
	return polygons

func fov():
	var fov_computer = get_node('FOV')
	var polygons = get_occluder_polygons()
	var fov_center = Vector2(0,0)
	var fov_poly = fov_computer.get_fov_from_polygons(polygons, fov_center)
	draw_colored_polygon(fov_poly, Color(0,0,255))
```

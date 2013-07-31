# Truss mesh SketchUp Plugin (2013)
#
# Saves lines as truss elements in VTK text format
# Author: Raul Durand 
# raul.durand@gmail.com
#
# Copy this file to SketchUp plugins path located at:
# [Program Files]\Google\Google SketchUp 8\Plugins\


# Function to erase all faces in model
def erase_faces()
	etys = Sketchup.active_model.entities
	faces = []
	for e in etys
		faces << e if e.is_a? Sketchup::Face
	end
	etys.erase_entities(faces)
end

# Draw point and edges ids
def draw_ids(points, edges)
	
	# Draw point ids
	for pt in points
		Sketchup.active_model.entities.add_text "n:"+pt.get_attribute("data","id").to_s, pt.position
	end

	# Draw edge ids
    edge_id = 0
	for edge in edges
		pt0 = edge.vertices[0].position
		pt1 = edge.vertices[1].position
		posx = 0.5*(pt0.x + pt1.x)
		posy = 0.5*(pt0.y + pt1.y)
		posz = 0.5*(pt0.z + pt1.z)
		pos = Geom::Point3d.new(posx, posy, posz)
		tag = edge.get_attribute("data", "tag", "")
		Sketchup.active_model.entities.add_text "e:"+edge_id.to_s+" "+tag, pos
		edge_id = edge_id + 1
	end

end

# Remove previous ids annotations
def remove_old_ids()
	etys = Sketchup.active_model.entities
	old_tags = []
	for ety in etys
		if ety.is_a? Sketchup::Text
			if ety.text[0..1]=="n:" or ety.text[0..1]=="e:"
				old_tags << ety
			end
		end 
	end

	# Remove old tags
	etys.erase_entities(old_tags)
end

# Set text tag for a selection
def set_text_tag()
	input = UI.inputbox ["Tag"], [""], "Tag for current selection"
	tag = input[0]
	sel = Sketchup.active_model.selection
	etys = Sketchup.active_model.entities

	# Set tags and print text
	old_tags = []
	for edge in sel
		if not edge.is_a? Sketchup::Edge
			next
		end

		# Set attribute
		edge.set_attribute("data", "tag", tag)

		# Insertion point
		pt0 = edge.vertices[0].position
		pt1 = edge.vertices[1].position
		posx = 0.5*(pt0.x + pt1.x)
		posy = 0.5*(pt0.y + pt1.y)
		posz = 0.5*(pt0.z + pt1.z)
		pos = Geom::Point3d.new(posx, posy, posz)

		# Remove tag text if any
		for ety in etys
			if ety.is_a? Sketchup::Text
				if pos == ety.point
					old_tags << ety
					break
				end
			end 
		end
		
		# Draw tag text
		Sketchup.active_model.entities.add_text "e:"+tag, pos

	end
	etys.erase_entities(old_tags)
end

def index_for_sort(point)
	x=point.x
	y=point.y
	z=point.z
	return 1000000*z + 1000*y + x
end

def edge_index_for_sort(edge)
	# Calc center of edge
	pt0 = edge.vertices[0].position
	pt1 = edge.vertices[1].position
	posx = 0.5*(pt0.x + pt1.x)
	posy = 0.5*(pt0.y + pt1.y)
	posz = 0.5*(pt0.z + pt1.z)
	return index_for_sort(Geom::Point3d.new(posx, posy, posz))
	
end

def get_geometric_etys()
	# Filter all edges in selection
	etys = Sketchup.active_model.entities.to_a
	edges = []
	for ety in etys
		if ety.is_a? Sketchup::Edge
			edges << ety
			if ety.get_attribute("data","tag")
			end
		end
	end

	edges.sort! {|a,b| edge_index_for_sort(a) <=> edge_index_for_sort(b) }	

	# Get vertices from edges
	points = []
	for edge in edges
		points = points | edge.vertices
	end

	# Join broken lines
	for point in points
		if point.edges.length>0
			ed1 = point.edges[0]
			ed2 = point.edges[0]
			# Incomplete code...
			# Needs to reconstruct edges array and reassign tag info
		end
	end

	# Sort vertices
	points.sort! {|a,b| index_for_sort(a.position) <=> index_for_sort(b.position) }

	# Set id as attribute 
	id = 0
	for point in points
		point.set_attribute "data", "id", id
		id = id+1
	end

	return points, edges

end

def draw_annotations()
	Sketchup.active_model.rendering_options['DisplayColorByLayer'] = true
	points, edges = get_geometric_etys()

	draw_ids points, edges

end


# Function to export lines as a truss mesh in VTK format
def export_as_vtk()
	
	# Change color by layer
	Sketchup.active_model.rendering_options['DisplayColorByLayer'] = true

	points, edges = get_geometric_etys()

	# Get filename from user interface
	filename = UI.savepanel "Enter the output file (vtk)", $default_path, $filename

	# Check for a valid filename
	begin 
		file = File.new filename, "w"
	rescue
		UI.messagebox('File could not be oppened for writting.')
		return
	end

	$default_path = File.dirname(filename)
	$filename     = File.basename(filename)

	units = Sketchup.active_model.options["UnitsOptions"]["LengthUnit"] # an index from 0 to 4
	baseunit = [1.inch, 1.feet, 1.mm, 1.cm, 1.m][units]
	factor = 1.0/baseunit

	# Print file header
	file.print "# vtk DataFile Version 3.0\n"     
	file.print "RbMesh output"  


	# Filter all edges in selection
	etys = Sketchup.active_model.entities.to_a
	has_tag = false
	for ety in etys
		if ety.is_a? Sketchup::Edge
			if ety.get_attribute("data","tag")
				has_tag = true
			end
		end
	end

	edges.sort! {|a,b| edge_index_for_sort(a) <=> edge_index_for_sort(b) }	

	# Get vertices from edges
	points = []
	for edge in edges
		points = points | edge.vertices
	end

	# Join broken lines
	for point in points
		if point.edges.length>0
			ed1 = point.edges[0]
			ed2 = point.edges[0]
			# Incomplete code...
			# Needs to reconstruct edges array and reassign tag info
		end
	end
	

	if has_tag
		# Find all tags
		tags = []
		for edge in edges
			tag = edge.get_attribute("data", "tag", "")
			tags = tags | [tag]
			edge.set_attribute("data", "tag_index", tags.index(tag))
		end

		# Print extra tag information	
		file.print " -- tags:"  
		print tags,"xxx\n"
		file.print tags[0]
		for tag in tags[1..-1]
			file.print ";" , tag
		end
	end

	file.print "\n"  
	file.print "ASCII\n" 
	file.print "DATASET UNSTRUCTURED_GRID\n"      
	file.print "\n"
	file.print "POINTS ", points.size,  " float\n" 

	# Print points
	for point in points
		file.printf "%15.3f %15.3f %15.3f \n", point.position.x*factor, point.position.y*factor, point.position.z*factor
	end
	file.printf "\n"

	# Write connectivities
	ndata = (2+1)*edges.size
	file.print "CELLS ", edges.size, " ", ndata, "\n"
	for edge in edges
		file.print 2, " "
		for point in edge.vertices
			file.print point.get_attribute("data", "id"), " " 
		end
		file.print "\n"
	end
	file.print "\n"

	# Write cell types
	file.print "CELL_TYPES ", edges.size, "\n"
	for cell in edges
		file.print "3 \n"
	end
	file.print "\n"

	# Write tags if available
	if has_tag

		file.print "CELL_DATA ", edges.size, "\n"
		file.print "SCALARS Tag int 1\n"
		file.print "LOOKUP_TABLE default\n"
		for cell in edges
			file.print cell.get_attribute("data","tag_index"),"\n"
		end
		file.print "\n"
	end

	file.close

	# Redraw ids
	remove_old_ids
	draw_ids points, edges

end


#UI.add_context_menu_handler {
#	|menu|
#	menu.add_separator
#	item = menu.add_item("Split") { split }
#	menu.set_validation_proc(item) { checkGroup }
#}


# Default path and filename
$default_path = "C:\\"
$filename     = "mesh.vtk"

# Plugin menu initialization
plugins_menu = UI.menu("Plugins")
item = plugins_menu.add_item("Set tag for selection") { set_text_tag }
item = plugins_menu.add_item("Draw annotations") { draw_annotations }
item = plugins_menu.add_item("Delete annotations") { remove_old_ids }
item = plugins_menu.add_item("Export mesh as .vtk") { export_as_vtk }


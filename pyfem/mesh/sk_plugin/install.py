
# Find Sketch Up Installation folder
import os

sk_path = os.environ['PROGRAMFILES'] + "\\" + "Google\\Google SketchUp 8\\Plugins" 
if not os.path.exists(sk_path):
    sk_path = os.environ['PROGRAMFILES(X86)'] + "\\" + "Google\\Google SketchUp 8\\Plugins" 

if not os.path.exists(sk_path):
    print "No SketchUp 8 installed"
    exit()

print sk_path
# Copying plugin
import shutil
shutil.copyfile('truss_plugin.rb', sk_path)

read_path="/projects/tallis/ekmolloy/multi-copy/data-old"
store_path="/projects/tallis/ekmolloy/multi-copy/data-old"

cd $read_path

read_files=$(ls -d */)
for f in $read_files
do
	cd $read_path/$f
	temp=$(ls -d */)
	for sub_file in $temp
	do
		cd $sub_file
		echo 1
		python /projects/tallis/ekmolloy/multi-copy/tools/generate_astral_map.py g_trees-raxml-sqlen-100.trees $store_path/$f/$sub_file/qikai_g_trees-raxml-sqlen-100_map.txt
		cd ..
	done
	cd ..
done

read_path="/projects/tallis/ekmolloy/multi-copy/data"
store_path="/projects/tallis/ekmolloy/multi-copy/qikaiy2-data"

cd $read_path

read_files=$(ls -d */)
for f in $read_files
do
	for i in {1..9}
	do
	python /projects/tallis/ekmolloy/multi-copy/tools/generate_astral_map.py $read_path/$f/0$i/g_trees-raxml-sqlen-100-ngen-100.trees $store_path/$f/0$i/g_trees-raxml-sqlen-100-ngen-100_map.txt
	python /projects/tallis/ekmolloy/multi-copy/tools/generate_astral_map.py $read_path/$f/0$i/g_trees-raxml-sqlen-100-ngen-25.trees  $store_path/$f/0$i/g_trees-raxml-sqlen-100-ngen-25_map.txt
	python /projects/tallis/ekmolloy/multi-copy/tools/generate_astral_map.py $read_path/$f/0$i/g_trees-raxml-sqlen-100-ngen-500.trees $store_path/$f/0$i/g_trees-raxml-sqlen-100-ngen-500_map.txt

	python /projects/tallis/ekmolloy/multi-copy/tools/generate_astral_map.py $read_path/$f/0$i/g_trees-raxml-sqlen-25-ngen-100.trees $store_path/$f/0$i/g_trees-raxml-sqlen-25-ngen-100_map.txt
	python /projects/tallis/ekmolloy/multi-copy/tools/generate_astral_map.py $read_path/$f/0$i/g_trees-raxml-sqlen-25-ngen-25.trees $store_path/$f/0$i/g_trees-raxml-sqlen-25-ngen-25_map.txt
	python /projects/tallis/ekmolloy/multi-copy/tools/generate_astral_map.py $read_path/$f/0$i/g_trees-raxml-sqlen-25-ngen-500.trees $store_path/$f/0$i/g_trees-raxml-sqlen-25-ngen-500_map.txt

	python /projects/tallis/ekmolloy/multi-copy/tools/generate_astral_map.py $read_path/$f/0$i/g_trees-raxml-sqlen-250-ngen-100.trees  $store_path/$f/0$i/g_trees-raxml-sqlen-250-ngen-100_map.txt
	python /projects/tallis/ekmolloy/multi-copy/tools/generate_astral_map.py $read_path/$f/0$i/g_trees-raxml-sqlen-250-ngen-25.trees $store_path/$f/0$i/g_trees-raxml-sqlen-250-ngen-25_map.txt
	python /projects/tallis/ekmolloy/multi-copy/tools/generate_astral_map.py $read_path/$f/0$i/g_trees-raxml-sqlen-250-ngen-500.trees $store_path/$f/0$i/g_trees-raxml-sqlen-250-ngen-500_map.txt

	python /projects/tallis/ekmolloy/multi-copy/tools/generate_astral_map.py $read_path/$f/0$i/g_trees-raxml-sqlen-50-ngen-100.trees $store_path/$f/0$i/g_trees-raxml-sqlen-50-ngen-100_map.txt
	python /projects/tallis/ekmolloy/multi-copy/tools/generate_astral_map.py $read_path/$f/0$i/g_trees-raxml-sqlen-50-ngen-25.trees $store_path/$f/0$i/g_trees-raxml-sqlen-50-ngen-25_map.txt
	python /projects/tallis/ekmolloy/multi-copy/tools/generate_astral_map.py $read_path/$f/0$i/g_trees-raxml-sqlen-50-ngen-500.trees $store_path/$f/0$i/g_trees-raxml-sqlen-50-ngen-500_map.txt
	done
	
done


import json
import argparse

def read_json(file_path):
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            data = json.load(file)
            return data
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except json.JSONDecodeError:
        print(f"Error: The file '{file_path}' is not a valid JSON file.")

def save_json(file_path, data):
    try:
        with open(file_path, 'w', encoding='utf-8') as file:
            json.dump(data, file, indent=4)
    except Exception as e:
        print(f"Error: Unable to save JSON file. {e}")
        

if "__main__" == __name__:
    p = argparse.ArgumentParser(description="Script to generate ACTS material map")
    p.add_argument(
        "-o",
        "--oldgeo",
        type=str,
        default="geometry-map.json",
        help="Output filename for the generated material map. Supported formats: JSON, CBOR.",
    )

    p.add_argument(
        "-n",
        "--newgeo",
        type=str,
        default="geometry-map-new.json",
        help="Output filename for the generated material map. Supported formats: JSON, CBOR.",
    )
    args = p.parse_args()

    file_path = args.oldgeo

    json_data = read_json(file_path)

    if json_data:
        surface = json_data["Surfaces"]["entries"]
        for auto in surface:
            if "approach" in auto and (auto["layer"] == 2 or auto["layer"] == 4 or auto["layer"] == 6 or auto["layer"] == 8 or auto["layer"] == 10):
                auto["value"]["material"]["binUtility"]["binningdata"][0]["bins"] = 20
                auto["value"]["material"]["binUtility"]["binningdata"][1]["bins"] = 40
                
                auto["value"]["material"]["mapMaterial"]=True

                """
                
                auto["value"]["material"]["mapMaterial"]=True
                if auto["approach"] == 1 and auto["layer"] != 10:
                    auto["value"]["material"]["binUtility"]["binningdata"][0]["bins"] = 314
                    auto["value"]["material"]["binUtility"]["binningdata"][1]["bins"] = 100
                    auto["value"]["material"]["mapMaterial"]=True
                elif auto["layer"] == 1 and auto["approach"] == 2:
                    auto["value"]["material"]["binUtility"]["binningdata"][0]["bins"] = 314
                    auto["value"]["material"]["binUtility"]["binningdata"][1]["bins"] = 100
                    auto["value"]["material"]["mapMaterial"]=True
                elif auto["layer"] == 18 and auto["approach"] == 2:
                    auto["value"]["material"]["binUtility"]["binningdata"][0]["bins"] = 314
                    auto["value"]["material"]["binUtility"]["binningdata"][1]["bins"] = 100
                    auto["value"]["material"]["mapMaterial"]=True
                elif auto["layer"] == 10 and auto["approach"] == 2:
                    auto["value"]["material"]["binUtility"]["binningdata"][0]["bins"] = 314
                    auto["value"]["material"]["binUtility"]["binningdata"][1]["bins"] = 100
                    auto["value"]["material"]["mapMaterial"]=True
                        """
                    
    save_json(args.newgeo, json_data)

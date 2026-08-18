import queue
from datetime import datetime
from time import sleep
import os


stdoutpipe = queue.Queue()
amigo = {"name": "AmiGo2", "collected": 0, "total": 0, "status": "not started", "last_update": "XX:XX"}
disgnet = {"name": "disgnet", "collected": 0, "total": 0, "status": "not started", "last_update": "XX:XX"}
hgnc = {"name": "HGNC", "collected": 0, "total": 0, "status": "not started", "last_update": "XX:XX"}
ensembl = {"name": "ensembl", "collected": 0, "total": 0, "status": "not started", "last_update": "XX:XX"}
gnomad = {"name": "gnomAD", "collected": 0, "total": 0, "status": "not started", "last_update": "XX:XX"}
vep = {"name": "VEP", "collected": 0, "total": 0, "status": "not started", "last_update": "XX:XX"}
running = True
output = queue.Queue()

def save_output(file_path):
    global output
    txt = ""

    while len(txt) < 10000:
        try:
            message = output.get_nowait()
        except queue.Empty:
            break

        txt += message
    if txt:
        with open(file_path, "a") as f:
            f.write(txt)

    return output.empty()

def send_message(mes, status=None, sender=None):
    global output, stdoutpipe
    time_stamp = datetime.now().strftime('%H%M')
    if sender is None and status is None:
        e = f"{time_stamp} Error\n{mes}\n"
        print(e)
        output.put(e)
    else:
        stdoutpipe.put((sender, status, mes, time_stamp))


def print_out():
    global amigo, disgnet, hgnc, ensembl, gnomad, vep
    time_stamp = datetime.now().strftime('%H%M')
    print(f"\n{time_stamp} Update Pipeline")
    for m in [amigo, disgnet, hgnc, ensembl, gnomad, vep]:
        print(f'{m["name"]} datasets: {m["collected"]}/{m["total"]} status: {m["last_update"]} {m["status"]}')


def sort_message(message):
    global output, running
    sender, status, information, ts = message
    match status:
        case 0:
            update_status(sender, information, ts)
            output.put(f"{ts} {sender} status update {information}\n")
        case 1:
            update_collected(sender, information, ts)
        case 2:
            update_total(sender, information, ts)
        case 3:
            e = f"{ts} results message: {information}\n"
            print(e)
            output.put(e)
        case _:
            e = f"{ts} {sender} unknown! status: {status} message: {information}\n"
            print(e)
            output.put(e)

def update_status(sender, information, ts):
    global amigo, disgnet, hgnc, ensembl, gnomad, vep
    match sender:
        case 1 | "Amigo" | "amigo" | "AmiGo" | "AmiGo2" | "amigo2" | "Amigo2":
            amigo["status"] = information
            amigo["last_update"] = ts
        case 2 | "disgnet":
            disgnet["status"] = information
            disgnet["last_update"] = ts
        case 3 | "hgnc" | "HGNC" | "hugo":
            hgnc["status"] = information
            hgnc["last_update"] = ts
        case 4 | "ensembl" | "ensemble" | "Ensembl" | "Ensemble":
            ensembl["status"] = information
            ensembl["last_update"] = ts
        case 5 | "gnomad" | "Gnomad" | "gnomAD" | "GnomAD":
            gnomad["status"] = information
            gnomad["last_update"] = ts
        case 6 | "VEP" | "vep":
            vep["status"] = information
            vep["last_update"] = ts
        case _ :
            e = f"{ts} {sender} unknown! message: {information}\n"
            print(e)
            output.put(e)


def update_total(sender, information, ts):
    global amigo, disgnet, hgnc, ensembl, gnomad, vep
    match sender:
        case 1 | "Amigo" | "amigo" | "AmiGo" | "AmiGo2" | "amigo2" | "Amigo2":
            amigo["total"] += information
            amigo["last_update"] = ts
        case 2 | "disgnet":
            disgnet["total"] += information
            disgnet["last_update"] = ts
        case 3 | "hgnc" | "HGNC" | "hugo":
            hgnc["total"] += information
            hgnc["last_update"] = ts
        case 4 | "ensembl":
            ensembl["total"] += information
            ensembl["last_update"] = ts
        case 5 | "gnomad" | "Gnomad" | "gnomAD" | "GnomAD":
            gnomad["total"] += information
            gnomad["last_update"] = ts
        case 6 | "VEP" | "vep":
            vep["total"] += information
            vep["last_update"] = ts
        case _ :
            e =f"{ts} {sender} unknown! message: {information}\n"
            print(e)
            output.put(e)

def update_collected(sender, information, ts):
    global amigo, disgnet, hgnc, ensembl, gnomad, vep
    match sender:
        case 1 | "Amigo" | "amigo" | "AmiGo" | "AmiGo2" | "amigo2" | "Amigo2":
            amigo["collected"] += information
            amigo["last_update"] = ts
        case 2 | "disgnet":
            disgnet["collected"] += information
            disgnet["last_update"] = ts
        case 3 | "hgnc" | "HGNC" | "hugo":
            hgnc["collected"] += information
            hgnc["last_update"] = ts
        case 4 | "ensembl":
            ensembl["collected"] += information
            ensembl["last_update"] = ts
        case 5 | "gnomad" | "Gnomad" | "gnomAD" | "GnomAD":
            gnomad["collected"] += information
            gnomad["last_update"] = ts
        case 6 | "VEP" | "vep":
            vep["collected"] += information
            vep["last_update"] = ts
        case _ :
            e = f"{ts} {sender} unknown! message: {information}\n"
            print(e)
            output.put(e)

def get_message():
    try:
        mes = stdoutpipe.get_nowait()
    except queue.Empty:
        return False

    sort_message(mes)
    return True

def stop():
    global running
    running = False
 
def run_io(file_path):
    if os.path.isfile(file_path):
        os.rename(file_path, f"{file_path}_old.txt")
    counter = 0
    run = True
    while run:
        if get_message():
            counter += 1
        else:
            sleep(10)
            counter += 100
        if counter > 1000:
            global running
            print_out()
            _ = save_output(file_path)
            counter = 0
            run = running

    while not save_output(file_path):
        print("saving output messages")


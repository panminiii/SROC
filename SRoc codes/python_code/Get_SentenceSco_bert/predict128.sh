cd /data/huangx
source ./myvenv/bin/activate
export BERT_BASE_DIR=/data/huangx/bert/bertbase/uncased_L-12_H-768_A-12
cd bert/bert
nohup python run_classifierPredict3.py --task_name=mnli --do_predict=true --data_dir=/data/huangx/bert/data128,/data/huangx/bert/dataout128 --vocab_file=$BERT_BASE_DIR/vocab.txt --bert_config_file=$BERT_BASE_DIR/bert_config.json --init_checkpoint=/data/huangx/bert/classfier128/bert_model.ckpt --max_seq_length=128 --output_dir=/data/huangx/bert/output128 &
